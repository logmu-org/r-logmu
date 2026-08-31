// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

// The R binding for veil. Everything R-specific lives here, never in src/veil/: this file owns the
// SEXPs and reads them out into plain C++ values before calling the core's builders. A Python or C#
// front end walks its own tree and calls the same builders, which is why the builder API -- not any
// wire format -- is the cross-language boundary.
//
// The R vocabulary ("+", "exp", "%%") is an R concern, so the name-to-Op mapping belongs here too;
// veil itself knows only its own monikers.

// MSVC C++ doesn't like the C99/C11 syntax R_ext/Complex.h uses for `Rcomplex`
#if defined(_MSC_VER) && !defined(R_LEGACY_RCOMPLEX)
#define R_LEGACY_RCOMPLEX
#endif

#include <cpp11.hpp>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>
#include "datey.h" // yearsFromClicks, clicksFromYears
#include "veil/AevRecipe.hpp"
#include "veil/Block.hpp"
#include "veil/ColumnScan.hpp"
#include "veil/ColumnSet.hpp"
#include "veil/ColumnView.hpp"
#include "veil/Engine.hpp"
#include "veil/Instruction.hpp"
#include "veil/Interpreter.hpp"
#include "veil/Node.hpp"
#include "veil/ObjStore.hpp"
#include "veil/Op.hpp"
#include "veil/RecordChunk.hpp"
#include "veil/StringMapping.hpp"
#include "veil/TimeGrid.hpp"
#include "veil/Tree.hpp"
#include "veil/TypeFull.hpp"

// The tree passes, in the order the pipeline runs them further down.
#include "veil/passAnnotateTypes.hpp"
#include "veil/passEncodeText.hpp"
#include "veil/passInsertCoercions.hpp"
#include "veil/passTagTimeVarying.hpp"
#include "veil/passCheckTimeInvariance.hpp"
#include "veil/passFoldConstants.hpp"
#include "veil/passNarrowComparisons.hpp"
#include "veil/passCheckSimilarityRange.hpp"
#include "veil/passPropagateIntervals.hpp"
#include "veil/passFoldIndicatorSquares.hpp"
#include "veil/passHoistFromIntegrate.hpp"
#include "veil/passShareCommonSubtrees.hpp"
#include "veil/passFoldIntervalComparisons.hpp"

// Lowering, which runs after all of those and ends the tree phase.
#include "veil/passLowerToBlock.hpp"

// The TAC passes, which run over the lowered block.
#include "veil/passAllocateVectorBuffers.hpp"
#include "veil/passComputeLiveness.hpp"
#include "veil/passComputeSlotDemand.hpp"
#include "veil/passFindDeadInstructions.hpp"

namespace
{

// Maps a canonical R head to a veil op. `it_ast` has already normalised the idiomatic aliases
// (`&&` to `&`, `pmax` to `max`), so only canonical heads arrive here.
//
// Arity disambiguates the two heads that carry both a unary and a binary meaning. `%in%` is absent
// because it desugars before its arguments are ingested -- see ingestIn.
std::optional<veil::Op> opFromRName(std::string_view name, size_t arity)
{
    if (name == "+") { return arity == 1 ? veil::Op::Pos : veil::Op::Add; }
    if (name == "-") { return arity == 1 ? veil::Op::Neg : veil::Op::Sub; }
    if (name == "*") { return veil::Op::Mul; }
    if (name == "/") { return veil::Op::RDiv; }
    if (name == "^") { return veil::Op::Pow; }
    if (name == "%%") { return veil::Op::IMod; }
    if (name == "%/%") { return veil::Op::IDiv; }

    if (name == "==") { return veil::Op::Eq; }
    if (name == "!=") { return veil::Op::Ne; }
    if (name == "<") { return veil::Op::Lt; }
    if (name == "<=") { return veil::Op::Le; }
    if (name == ">") { return veil::Op::Gt; }
    if (name == ">=") { return veil::Op::Ge; }

    if (name == "&") { return veil::Op::And; }
    if (name == "|") { return veil::Op::Or; }
    if (name == "!") { return veil::Op::Not; }
    if (name == "xor") { return veil::Op::Xor; }
    if (name == "is.na") { return veil::Op::IsNa; }

    if (name == "abs") { return veil::Op::Abs; }
    if (name == "sqrt") { return veil::Op::Sqrt; }
    if (name == "exp") { return veil::Op::Exp; }
    if (name == "log") { return veil::Op::Log; }
    if (name == "log10") { return veil::Op::Log10; }
    if (name == "log1p") { return veil::Op::Log1p; }
    if (name == "expm1") { return veil::Op::Expm1; }
    if (name == "sin") { return veil::Op::Sin; }
    if (name == "cos") { return veil::Op::Cos; }
    if (name == "floor") { return veil::Op::Floor; }
    if (name == "ceiling") { return veil::Op::Ceiling; }
    if (name == "round") { return veil::Op::Round; }
    if (name == "trunc") { return veil::Op::Trunc; }
    if (name == "sign") { return veil::Op::Sign; }

    if (name == "min") { return veil::Op::Min; }
    if (name == "max") { return veil::Op::Max; }
    if (name == "clamp") { return veil::Op::Clamp; }

    // Both spellings of the same thing. An `if` without an else arrives with arity 2 and is
    // rejected by the arity check: a value is wanted on both paths.
    if (name == "if" || name == "ifelse") { return veil::Op::Select; }

    return std::nullopt;
}

veil::ColumnId resolveColumn(const cpp11::strings& columnNames, const std::string& name)
{
    for (R_xlen_t i = 0; i < columnNames.size(); ++i)
    {
        if (static_cast<std::string>(columnNames[i]) == name)
        {
            return static_cast<veil::ColumnId>(i);
        }
    }
    cpp11::stop("Field `.i$%s` is not a column of the data.", name.c_str());
}

// A folded constant crosses as whatever R type it folded to.
//
// datey and durationy are both 32-bit integer clicks on the same annual grid, so an integer literal
// carries the right number either way and converts to years by the same division. Only the type tag
// differs, and it is the R class that distinguishes them -- which matters once the type rules bite,
// since `datey - datey` is a durationy and `datey + datey` is nonsense.
//
// A bare integer with no class is a count, and becomes a double: raw integer is not a veil type
// (integer columns are rejected precisely because count-versus-category is ambiguous).
veil::NodeId ingestLit(veil::Tree& tree, SEXP value)
{
    if (Rf_xlength(value) != 1)
    {
        cpp11::stop(
            "A veil literal must be a single value, not a vector of length %lld.",
            static_cast<long long>(Rf_xlength(value)));
    }

    switch (TYPEOF(value))
    {
        case LGLSXP: return tree.buildLitBool(cpp11::as_cpp<bool>(value));
        case REALSXP: return tree.buildLitDouble(cpp11::as_cpp<double>(value));
        case STRSXP: return tree.buildLitText(cpp11::as_cpp<std::string>(value));
        case INTSXP:
        {
            const int clicks = cpp11::as_cpp<int>(value);
            if (Rf_inherits(value, "durationy")) { return tree.buildLitInt(clicks, veil::TypeFull::createDurationy()); }
            if (Rf_inherits(value, "datey")) { return tree.buildLitInt(clicks, veil::TypeFull::createDatey()); }
            return tree.buildLitDouble(static_cast<double>(clicks));
        }
    }

    cpp11::stop("A veil literal of R type '%s' is not readable.", Rf_type2char(TYPEOF(value)));
}

veil::NodeId ingest(
    veil::Tree& tree,
    veil::ObjStore& objs,
    const cpp11::list& node,
    const cpp11::strings& columnNames);

int scalarClicks(SEXP value, const char* what)
{
    if (Rf_xlength(value) != 1 || TYPEOF(value) != INTSXP)
    {
        cpp11::stop("A veil %s bound must be a single datey/durationy value.", what);
    }
    const int clicks = cpp11::as_cpp<int>(value);
    if (clicks == NA_INTEGER)
    {
        // Refused at the crossing rather than left to mean INT_MIN further in, where it would become
        // a table origin or an interval bound four millennia adrift.
        cpp11::stop("A veil %s bound is missing.", what);
    }
    return clicks;
}

// Reads a mortality_table out of R into the core's own representation. The log mu matrix is copied
// verbatim, so it keeps R's column-major order -- which is what MortalityTable documents and
// indexes by.
veil::ObjId readTable(veil::ObjStore& objs, SEXP obj)
{
    const cpp11::doubles logMu = cpp11::as_cpp<cpp11::doubles>(obj);
    const cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(Rf_getAttrib(obj, R_DimSymbol));
    if (dim.size() != 2)
    {
        cpp11::stop("A mortality_table must hold a 2D age-period matrix.");
    }

    veil::MortalityTable table;
    table.ageCount = static_cast<uint32_t>(dim[0]);
    table.periodCount = static_cast<uint32_t>(dim[1]);
    table.x0Clicks = scalarClicks(Rf_getAttrib(obj, Rf_install("x0")), "x0");
    table.t0Clicks = scalarClicks(Rf_getAttrib(obj, Rf_install("t0")), "t0");
    table.logMu.assign(logMu.begin(), logMu.end());

    if (!table.isWellFormed())
    {
        cpp11::stop("A mortality_table's log_mu does not match its dimensions.");
    }

    return objs.add(std::move(table));
}

veil::ObjId readInclude(
    veil::Tree& tree,
    veil::ObjStore& objs,
    SEXP obj,
    const cpp11::strings& columnNames)
{
    veil::Include include;

    // AN INDICATOR IS AN INCLUDE THAT CARRIES NO `terms`. It holds an `ast` and means "this gate
    // must hold, over all of time", which is one gate term and nothing else.
    //
    // THIS IS THE ORDINARY SHAPE, NOT A CORNER: `include(.i$sex == "male")` returns an indicator by
    // design, so it is what a user's population include arrives as. Handled here rather than by
    // having R rewrite it first, because the caller's guard promises to accept any `include` object
    // and an indicator is one -- leaving it out failed with "expected 'list' actual 'NULL'", which
    // names neither the argument nor the cause.
    if (Rf_inherits(obj, "indicator"))
    {
        const cpp11::list indicator = cpp11::as_cpp<cpp11::list>(obj);
        const veil::NodeId gate =
            ingest(tree, objs, cpp11::as_cpp<cpp11::list>(indicator["ast"]), columnNames);
        include.terms.push_back(veil::GateTerm{gate});
        return objs.add(std::move(include));
    }

    const cpp11::list includeList = cpp11::as_cpp<cpp11::list>(obj);
    const cpp11::list terms = cpp11::as_cpp<cpp11::list>(includeList["terms"]);

    include.terms.reserve(static_cast<size_t>(terms.size()));

    for (R_xlen_t i = 0; i < terms.size(); ++i)
    {
        const cpp11::list term = cpp11::as_cpp<cpp11::list>(terms[i]);
        const std::string kind = cpp11::as_cpp<std::string>(term["kind"]);

        if (kind == "absolute")
        {
            include.terms.push_back(veil::AbsoluteTerm{scalarClicks(term["from"], "absolute from"),
                                                       scalarClicks(term["to"], "absolute to")});
        }
        else if (kind == "offset")
        {
            // Ingested into this same tree, exactly as a gate is: the offset is an expression, and
            // the passes that matter for it sweep every node rather than walking down from a root.
            const veil::NodeId offset =
                ingest(tree, objs, cpp11::as_cpp<cpp11::list>(term["offset"]), columnNames);
            include.terms.push_back(veil::OffsetTerm{offset,
                                                     scalarClicks(term["from"], "offset from"),
                                                     scalarClicks(term["to"], "offset to")});
        }
        else if (kind == "gate")
        {
            const veil::NodeId gate = ingest(tree, objs, cpp11::as_cpp<cpp11::list>(term["ast"]), columnNames);
            include.terms.push_back(veil::GateTerm{gate});
        }
        else
        {
            cpp11::stop("Unknown include term kind '%s'.", kind.c_str());
        }
    }

    return objs.add(std::move(include));
}

// Lowers a spliced concept object to the parameters the core needs, per the object-leaf rules:
// a const mortality becomes a scalar; a table becomes `vector_log_mu(table, birth)`; a mortality
// expression splices its own tree; an include crosses as bounds plus gates.
veil::NodeId ingestObj(
    veil::Tree& tree,
    veil::ObjStore& objs,
    SEXP obj,
    const cpp11::strings& columnNames)
{
    if (Rf_inherits(obj, "mortality_const"))
    {
        return tree.buildLitDouble(cpp11::as_cpp<double>(obj));
    }

    if (Rf_inherits(obj, "mortality_table"))
    {
        const veil::ObjId table = readTable(objs, obj);
        const veil::NodeId tableNode = tree.buildObj(table);
        const veil::NodeId birth = tree.buildField(resolveColumn(columnNames, "birth"), "birth");
        return tree.buildCall(veil::Op::VectorLogMu, {tableNode, birth});
    }

    if (Rf_inherits(obj, "mortality_expr"))
    {
        const cpp11::list expr = cpp11::as_cpp<cpp11::list>(obj);
        return ingest(tree, objs, cpp11::as_cpp<cpp11::list>(expr["ast"]), columnNames);
    }

    // A VARIABLE SPLICES ITS OWN TREE, exactly as a mortality expression above does: both are a
    // `list(ast = ...)` behind a class, and both mean "evaluate this expression here".
    //
    // REACHED ONLY BY A DIRECT OBJ LEAF, not through `aev()`. R's `it_classify_value` splices the
    // `ast` of any concept object that has one, so a `variable` written into a pronoun expression
    // never arrives here -- its tree does. This exists so that the entry point accepts an obj leaf
    // of any concept type rather than only the ones R happens not to splice, and it is pinned by a
    // test that hands it one directly. The same is true of the `mortality_expr` branch above.
    //
    // Placed before `include` because an `indicator` is BOTH a variable and an include: as a leaf in
    // a value position it is being read for its {0,1} value. No test distinguishes the two orders,
    // because a single-gate include lowers back to that gate's value -- so this is the reading that
    // is right rather than the one that is observable.
    if (Rf_inherits(obj, "variable"))
    {
        const cpp11::list expr = cpp11::as_cpp<cpp11::list>(obj);
        return ingest(tree, objs, cpp11::as_cpp<cpp11::list>(expr["ast"]), columnNames);
    }

    if (Rf_inherits(obj, "include"))
    {
        return tree.buildObj(readInclude(tree, objs, obj, columnNames));
    }

    cpp11::stop("A veil object leaf of this class cannot be lowered.");
}

// One element of an `%in%` set. The whitelist admits character and numeric sets only; an unclassed
// integer is a count, as it is anywhere else.
veil::NodeId buildSetElement(veil::Tree& tree, SEXP set, R_xlen_t i)
{
    switch (TYPEOF(set))
    {
        // STRING_ELT yields a CHARSXP, which as_cpp<std::string> will not take -- it wants a
        // length-1 STRSXP. r_string is the conversion that handles a single element.
        case STRSXP: return tree.buildLitText(static_cast<std::string>(cpp11::r_string(STRING_ELT(set, i))));
        case REALSXP: return tree.buildLitDouble(REAL(set)[i]);
        case INTSXP: return tree.buildLitDouble(static_cast<double>(INTEGER(set)[i]));
        case LGLSXP: return tree.buildLitBool(LOGICAL(set)[i] != 0);
    }
    cpp11::stop(
        "The right-hand side of `%%in%%` must be a character or numeric set, not '%s'.",
        Rf_type2char(TYPEOF(set)));
}

// `%in%` is the one place the whitelist admits a vector literal, and it means what writing the
// disjunction by hand would mean: `x %in% c(a, b)` is `x == a | x == b`. Desugaring it here keeps
// the op set unchanged and needs no set-valued operand.
//
// It has to run before the arguments are ingested, because the right-hand side is a vector and
// ingestLit rightly refuses anything that is not a single value.
//
// The left-hand side is ingested once and shared by every comparison, so the result is a DAG rather
// than a tree. That is the shape CSE would produce anyway, and nothing downstream assumes otherwise.
//
// A text set desugars the same way and costs nothing extra: passEncodeText resolves each element on
// its own, so an element naming no level settles its own comparison and drops out. Turning the whole
// chain into one integer-set test held in the ObjStore would be O(1) rather than O(n) in the set
// size, which at two or three categories is worth nothing; the tree phase can recognise the pattern
// later if it ever is.
veil::NodeId ingestIn(
    veil::Tree& tree,
    veil::ObjStore& objs,
    const cpp11::list& node,
    const cpp11::strings& columnNames)
{
    const cpp11::list args = cpp11::as_cpp<cpp11::list>(node["args"]);
    if (args.size() != 2)
    {
        cpp11::stop("`%%in%%` takes two arguments, but %d were given.", static_cast<int>(args.size()));
    }

    const veil::NodeId lhs = ingest(tree, objs, cpp11::as_cpp<cpp11::list>(args[0]), columnNames);

    const cpp11::list rhs = cpp11::as_cpp<cpp11::list>(args[1]);
    if (cpp11::as_cpp<std::string>(rhs["kind"]) != "lit")
    {
        cpp11::stop("The right-hand side of `%%in%%` must be a literal set.");
    }

    SEXP set = rhs["value"];
    const R_xlen_t count = Rf_xlength(set);
    if (count == 0)
    {
        return tree.buildLitBool(false); // Membership of nothing is false, never NA.
    }

    veil::NodeId acc = veil::invalidNodeId;
    for (R_xlen_t i = 0; i < count; ++i)
    {
        const veil::NodeId eq = tree.buildCall(veil::Op::Eq, {lhs, buildSetElement(tree, set, i)});
        acc = (i == 0) ? eq : tree.buildCall(veil::Op::Or, {acc, eq});
    }
    return acc;
}

// Walks the it_node tree, building the core's own nodes as it goes. A call builds its children
// first, then itself, so a child's NodeId always exists before its parent refers to it.
veil::NodeId ingest(
    veil::Tree& tree,
    veil::ObjStore& objs,
    const cpp11::list& node,
    const cpp11::strings& columnNames)
{
    const std::string kind = cpp11::as_cpp<std::string>(node["kind"]);

    if (kind == "time") { return tree.buildTime(); }
    if (kind == "lit") { return ingestLit(tree, node["value"]); }

    if (kind == "field")
    {
        const std::string name = cpp11::as_cpp<std::string>(node["name"]);
        return tree.buildField(resolveColumn(columnNames, name), name);
    }

    if (kind == "call")
    {
        // `%in%` is desugared, so it is intercepted before its arguments are ingested.
        if (cpp11::as_cpp<std::string>(node["fn"]) == "%in%")
        {
            return ingestIn(tree, objs, node, columnNames);
        }

        const cpp11::list args = cpp11::as_cpp<cpp11::list>(node["args"]);
        std::vector<veil::NodeId> childIds;
        childIds.reserve(static_cast<size_t>(args.size()));
        for (R_xlen_t i = 0; i < args.size(); ++i)
        {
            childIds.push_back(ingest(tree, objs, cpp11::as_cpp<cpp11::list>(args[i]), columnNames));
        }

        const std::string fn = cpp11::as_cpp<std::string>(node["fn"]);
        const std::optional<veil::Op> op = opFromRName(fn, childIds.size());
        if (!op.has_value())
        {
            cpp11::stop("`%s` has no veil operation.", fn.c_str());
        }

        const veil::OpInfo& info = veil::opInfo(*op);
        if (!info.implemented)
        {
            cpp11::stop(
                "The veil operation `%s` is declared but not yet implemented.",
                std::string(info.moniker).c_str());
        }
        if (info.arity != childIds.size())
        {
            cpp11::stop("`%s` takes %d argument(s), but %d were given.",
                        fn.c_str(), static_cast<int>(info.arity), static_cast<int>(childIds.size()));
        }

        return tree.buildCall(*op, std::move(childIds));
    }

    if (kind == "obj")
    {
        return ingestObj(tree, objs, node["value"], columnNames);
    }

    cpp11::stop("Unknown it_node kind '%s'.", kind.c_str());
}

// The crossing's strings, and what each factor column's own level codes become.
//
// R INTERNS, VEIL MERGES. Turning a million strings into codes is what `factor()` already does, in C,
// so the binding never touches a row's text; what it does is put the level tables of every text
// column onto ONE numbering, which is a few strings per column and costs nothing. That is the whole
// division of labour, and it is why a literal means the same thing in every column and every dataset
// of a batch without anything having to be reconciled later.
struct TextEncoding final
{
    veil::StringMapping mapping;

    // One entry per column, empty unless that column is a factor -- when it maps the factor's own
    // level numbering onto the shared one, indexed from zero as R's codes are indexed from one.
    std::vector<std::vector<int32_t>> levelIndices;

    int levels() const { return static_cast<int>(this->mapping.size()); }

    const std::vector<int32_t>& indicesFor(size_t column) const { return this->levelIndices[column]; }
};

TextEncoding buildTextEncoding(cpp11::list columns, const cpp11::strings& columnNames)
{
    TextEncoding encoding;
    encoding.levelIndices.resize(static_cast<size_t>(columns.size()));

    for (R_xlen_t i = 0; i < columns.size(); ++i)
    {
        SEXP column = VECTOR_ELT(SEXP(columns), i);
        if (TYPEOF(column) != INTSXP || !Rf_inherits(column, "factor")) { continue; }

        SEXP levels = Rf_getAttrib(column, R_LevelsSymbol);
        if (TYPEOF(levels) != STRSXP)
        {
            cpp11::stop("Column `%s` is a factor with no levels attribute.",
                        static_cast<std::string>(columnNames[i]).c_str());
        }

        std::vector<int32_t>& indices = encoding.levelIndices[static_cast<size_t>(i)];
        indices.reserve(static_cast<size_t>(Rf_xlength(levels)));
        for (R_xlen_t level = 0; level < Rf_xlength(levels); ++level)
        {
            indices.push_back(
                encoding.mapping.addString(static_cast<std::string>(cpp11::r_string(STRING_ELT(levels, level)))));
        }
    }

    return encoding;
}

// The veil type of an R column, mirroring the readability rules and ingestLit's dispatch. The classed
// integer cases are checked before the bare-integer rejection, since a factor and a datey are both
// INTSXP with a class.
//
// A FACTOR'S TYPE CARRIES THE CROSSING'S LEVEL COUNT, NOT ITS OWN. Its codes have been remapped onto
// the shared numbering, so its own count would understate the range and the scan would read a valid
// code as a missing one.
veil::TypeFull rColumnType(SEXP column, const std::string& name, int categoryLevels)
{
    switch (TYPEOF(column))
    {
        case LGLSXP: return veil::TypeFull::createBool();
        case REALSXP: return veil::TypeFull::createDouble();
        case STRSXP: return veil::TypeFull::createText();
        case INTSXP:
        {
            if (Rf_inherits(column, "durationy")) { return veil::TypeFull::createDurationy(); }
            if (Rf_inherits(column, "datey")) { return veil::TypeFull::createDatey(); }
            if (Rf_inherits(column, "factor")) { return veil::TypeFull::createCategory(categoryLevels); }
            cpp11::stop("Column `%s` is a bare integer; integer columns are not readable "
                        "(count versus category is ambiguous).", name.c_str());
        }
    }
    cpp11::stop("Column `%s` has type '%s', which veil cannot read.", name.c_str(),
                Rf_type2char(TYPEOF(column)));
}

// The types the reader and scan can handle so far. Text is absent because it is not a column type:
// a column of strings is a factor, and arrives as a category.
bool isScannable(veil::Type type)
{
    return type == veil::Type::Double || type == veil::Type::Datey
        || type == veil::Type::Durationy || type == veil::Type::Bool
        || type == veil::Type::Category;
}

cpp11::strings columnNamesOf(cpp11::list columns)
{
    const SEXP namesSexp = Rf_getAttrib(SEXP(columns), R_NamesSymbol);
    if (namesSexp == R_NilValue)
    {
        cpp11::stop("`columns` must be a named list of data columns.");
    }
    return cpp11::strings(namesSexp);
}

// Materialises one R column into a ColumnSet: a zero-copy view for the types that need no copy, a
// converted buffer for a logical. The view points into R's memory, so `column` must stay protected
// for as long as the set is read. Returns the id of the added column.
veil::ColumnId readColumn(
    veil::ColumnSet& set,
    const std::string& name,
    SEXP column,
    veil::TypeFull type,
    const std::vector<int32_t>& levelIndices)
{
    const size_t count = static_cast<size_t>(Rf_xlength(column));
    switch (type.type)
    {
        case veil::Type::Double:
            return set.addView(name, type, REAL(column), count);
        case veil::Type::Datey:
        case veil::Type::Durationy:
            return set.addView(name, type, INTEGER(column), count);
        case veil::Type::Bool:
        {
            // R holds a logical as int32; convert to a 1-byte bool the set owns. A missing logical is
            // refused, since a logical column must not contain NA.
            auto values = std::make_unique<bool[]>(count);
            const int* raw = LOGICAL(column);
            for (size_t i = 0; i < count; ++i)
            {
                if (raw[i] == NA_LOGICAL)
                {
                    cpp11::stop(
                        "Column `%s` has a missing logical value; a logical column must not "
                        "contain NA.",
                        name.c_str());
                }
                values[i] = raw[i] != 0;
            }
            return set.addBools(name, std::move(values), count);
        }
        case veil::Type::Category:
        {
            // R numbers a factor's levels from one; the core numbers the crossing's strings from
            // zero. This copy puts every column's codes on that one shared numbering, so a code means
            // the same thing whichever column it came from and whichever front end supplied it.
            //
            // A MISSING LEVEL IS REFUSED, exactly as a missing logical is a few lines above. A
            // category is stored as an `int`, and an `int` has no platform-independent NA
            // representation to carry absence in -- R's NA_INTEGER is INT_MIN, which no other
            // platform repeats, and a sentinel invented here would be no more portable for being
            // ours. Only `double`, `datey` and `durationy` have a missing state, and they have it
            // because their own representations provide one.
            //
            // An anything-but-a-level code is caught by the same test, so the remap cannot read
            // outside `levelIndices` and every code handed on is a genuine index into the mapping.
            auto values = std::make_unique<int32_t[]>(count);
            const int* raw = INTEGER(column);
            for (size_t i = 0; i < count; ++i)
            {
                const int level = raw[i];
                if (level < 1 || static_cast<size_t>(level) > levelIndices.size())
                {
                    cpp11::stop(
                        "Column `%s` has a missing value; a factor column must not contain NA. Give "
                        "the unknown values a level of their own, with `addNA()` or an explicit "
                        "level such as \"unknown\".",
                        name.c_str());
                }
                values[i] = levelIndices[static_cast<size_t>(level) - 1];
            }
            return set.addCodes(name, std::move(values), count, type.max);
        }
        default:
            cpp11::stop("Column `%s`: veil cannot read a column of this type. A column of strings "
                        "must be a factor.", name.c_str());
    }
}

// Reads every column into `set` and records each one's type and, where it can be scanned, its
// constraints. A column the scan cannot read is left unscanned rather than refused, so it does not
// stop the rest of an expression being specialised. `set` is filled rather than returned because the
// views inside it point at memory it owns, so it must outlive every use.
//
// Defined after readColumn, which it calls: this file has no forward declarations for its helpers.
// Returns the crossing's TextEncoding, which it has to build first in any case: a factor's veil type
// depends on the shared numbering, so the merge happens before the first column is typed. The caller
// keeps it because the mapping outlives this call -- passEncodeText resolves literals against it.
TextEncoding prepareColumns(
    cpp11::list columns,
    const cpp11::strings& columnNames,
    veil::ColumnSet& set,
    std::vector<veil::TypeFull>& columnTypes,
    std::vector<std::optional<veil::TypeWithConstraints>>& constraints)
{
    TextEncoding encoding = buildTextEncoding(columns, columnNames);

    columnTypes.reserve(static_cast<size_t>(columns.size()));
    constraints.reserve(static_cast<size_t>(columns.size()));

    for (R_xlen_t i = 0; i < columns.size(); ++i)
    {
        const std::string name = static_cast<std::string>(columnNames[i]);
        SEXP column = VECTOR_ELT(SEXP(columns), i);
        const veil::TypeFull type = rColumnType(column, name, encoding.levels());
        columnTypes.push_back(type);
        if (isScannable(type.type))
        {
            const veil::ColumnId id =
                readColumn(set, name, column, type, encoding.indicesFor(static_cast<size_t>(i)));
            constraints.push_back(veil::scanColumn(set.at(id)));
        }
        else
        {
            constraints.push_back(std::nullopt);
        }
    }

    return encoding;
}

} // namespace

const char* typeName(const veil::TypeFull& type)
{
    switch (type.type)
    {
        case veil::Type::Bool: return "bool";
        case veil::Type::Double: return "double";
        case veil::Type::Datey: return "datey";
        case veil::Type::DateyInterval: return "datey_interval";
        case veil::Type::Durationy: return "durationy";
        case veil::Type::Text: return "text";
        case veil::Type::Category: return "category";
    }
    return "unknown";
}

using namespace cpp11::literals;

cpp11::list dumpObj(const veil::Obj& obj)
{
    if (const auto* c = std::get_if<veil::MortalityConst>(&obj))
    {
        return cpp11::writable::list({"kind"_nm = "mortality_const", "log_mu"_nm = c->logMu});
    }

    if (const auto* t = std::get_if<veil::MortalityTable>(&obj))
    {
        // log_mu goes back in the order the core holds it, so a test can compare it against the
        // source matrix directly and catch a column-major mix-up.
        cpp11::writable::doubles logMu(static_cast<R_xlen_t>(t->logMu.size()));
        for (size_t i = 0; i < t->logMu.size(); ++i) { logMu[static_cast<R_xlen_t>(i)] = t->logMu[i]; }

        return cpp11::writable::list({
            "kind"_nm = "mortality_table",
            "age_count"_nm = static_cast<int>(t->ageCount),
            "period_count"_nm = static_cast<int>(t->periodCount),
            "x0_clicks"_nm = t->x0Clicks,
            "t0_clicks"_nm = t->t0Clicks,
            "log_mu"_nm = logMu,
        });
    }

    const auto& inc = std::get<veil::Include>(obj);
    cpp11::writable::list terms;
    for (const auto& term : inc.terms)
    {
        if (const auto* a = std::get_if<veil::AbsoluteTerm>(&term))
        {
            terms.push_back(cpp11::writable::list({
                "kind"_nm = "absolute", "from_clicks"_nm = a->fromClicks, "to_clicks"_nm = a->toClicks}));
        }
        else if (const auto* o = std::get_if<veil::OffsetTerm>(&term))
        {
            terms.push_back(cpp11::writable::list({
                "kind"_nm = "offset", "offset"_nm = static_cast<int>(o->offset),
                "from_clicks"_nm = o->fromClicks, "to_clicks"_nm = o->toClicks}));
        }
        else
        {
            terms.push_back(cpp11::writable::list({
                "kind"_nm = "gate", "ast"_nm = static_cast<int>(std::get<veil::GateTerm>(term).ast)}));
        }
    }

    return cpp11::writable::list({"kind"_nm = "include", "terms"_nm = terms});
}

// A literal's value rendered for inspection from R. Folding passes decide values, not just shapes, so
// a dump that showed only kinds and types could not tell a fold to `true` from a fold to `false`.
// A CATEGORY LITERAL RENDERS AS ITS STRING, NOT AS ITS INDEX, which is what the crossing's
// StringMapping is kept around to do. `.i$sex == "male"` becomes an integer comparison before
// anything reads this, so without the mapping the only account of it anyone could get back would be
// against a `0` whose meaning lives in a table they cannot see.
//
// This is the only direction the core reads the mapping in, and it reads it only while explaining
// itself. Nothing on a calculating path renders a string.
std::string formatLiteral(const veil::LitPayload& lit, const veil::StringMapping& mapping)
{
    if (const auto* b = std::get_if<bool>(&lit.value)) { return *b ? "TRUE" : "FALSE"; }
    if (const auto* i = std::get_if<int>(&lit.value))
    {
        if (lit.type.type == veil::Type::Category) { return mapping.stringAtIndex(*i); }
        return std::to_string(*i);
    }
    if (const auto* s = std::get_if<std::string>(&lit.value)) { return *s; }

    const double value = std::get<double>(lit.value);
    std::ostringstream out;
    out << std::setprecision(17) << value; // Enough to round-trip, without a fixed-point tail of zeros.
    return out.str();
}

// The per-node dump shared by the tree-pipeline entry points: each node's kind, resolved type, a
// label (a moniker for a call, a name for a field), the value of a literal, and whether it is
// time-varying. `intervals` is optional; when no pass has computed them the bounds report NA, so
// every entry point returns the same shape whether or not it ran interval propagation.
cpp11::list dumpTree(
    const veil::Tree& tree,
    const std::vector<char>& timeVarying,
    const veil::StringMapping& mapping,
    const std::vector<veil::Interval>* intervals = nullptr)
{
    cpp11::writable::strings kinds;
    cpp11::writable::strings types;
    cpp11::writable::strings labels;
    cpp11::writable::strings values;
    cpp11::writable::logicals vectors;
    cpp11::writable::doubles los;
    cpp11::writable::doubles his;
    for (veil::NodeId id = 0; id < static_cast<veil::NodeId>(tree.size()); ++id)
    {
        const veil::Node& n = tree.at(id);

        if (intervals != nullptr && id < intervals->size())
        {
            los.push_back((*intervals)[id].lo);
            his.push_back((*intervals)[id].hi);
        }
        else
        {
            los.push_back(NA_REAL);
            his.push_back(NA_REAL);
        }

        // A mortality table obj carries no value type, and reports an empty string.
        types.push_back(n.type.has_value() ? typeName(*n.type) : "");
        vectors.push_back(static_cast<cpp11::r_bool>(id < timeVarying.size() && timeVarying[id] != 0));
        values.push_back(veil::isLit(n) ? formatLiteral(std::get<veil::LitPayload>(n.payload), mapping) : "");

        if (veil::isCall(n))
        {
            kinds.push_back("call");
            labels.push_back(std::string(veil::moniker(std::get<veil::CallPayload>(n.payload).op)));
        }
        else if (veil::isField(n))
        {
            kinds.push_back("field");
            labels.push_back(std::get<veil::FieldPayload>(n.payload).name);
        }
        else
        {
            kinds.push_back(veil::isLit(n) ? "lit" : veil::isTime(n) ? "time" : "obj");
            labels.push_back("");
        }
    }

    return cpp11::writable::list({
        "node_count"_nm = static_cast<int>(tree.size()),
        "root"_nm = static_cast<int>(tree.root()),
        "kinds"_nm = kinds,
        "types"_nm = types,
        "labels"_nm = labels,
        "values"_nm = values,
        "vectors"_nm = vectors,
        "lo"_nm = los,
        "hi"_nm = his,
    });
}

// Builds a veil tree from an it_node and reports what was built -- the tree and the objects it
// refers to -- so the crossing can be tested from R before any of the passes downstream of it exist.
[[cpp11::register]]
cpp11::list cpp_veil_ingest_ast(cpp11::list node, cpp11::strings column_names)
{
    veil::Tree tree;
    veil::ObjStore objs;
    const veil::NodeId root = ingest(tree, objs, node, column_names);
    tree.setRoot(root);

    cpp11::writable::strings monikers;
    cpp11::writable::strings litTypes;
    cpp11::writable::strings kinds;
    for (veil::NodeId id = 0; id < static_cast<veil::NodeId>(tree.size()); ++id)
    {
        const veil::Node& n = tree.at(id);
        if (veil::isCall(n))
        {
            kinds.push_back("call");
            monikers.push_back(std::string(veil::moniker(std::get<veil::CallPayload>(n.payload).op)));
        }
        else if (veil::isLit(n))
        {
            kinds.push_back("lit");
            litTypes.push_back(typeName(std::get<veil::LitPayload>(n.payload).type));
        }
        else if (veil::isField(n)) { kinds.push_back("field"); }
        else if (veil::isTime(n)) { kinds.push_back("time"); }
        else { kinds.push_back("obj"); }
    }

    cpp11::writable::list objDump;
    for (veil::ObjId id = 0; id < static_cast<veil::ObjId>(objs.size()); ++id)
    {
        objDump.push_back(dumpObj(objs.at(id)));
    }

    return cpp11::writable::list({
        "node_count"_nm = static_cast<int>(tree.size()),
        "root"_nm = static_cast<int>(tree.root()),
        "kinds"_nm = kinds,
        "monikers"_nm = monikers,
        "lit_types"_nm = litTypes,
        "objs"_nm = objDump,
    });
}

// Builds a veil tree, annotates every node's type, and reports the per-node kind, type and label, so
// the type annotation pass can be tested from R. `columns` is a named list of the actual data
// columns; only their types are read here, though passing the data now keeps the signature stable
// for the constraint scan that reads values later.
[[cpp11::register]]
cpp11::list cpp_veil_prepare(cpp11::list node, cpp11::list columns)
{
    const cpp11::strings columnNames = columnNamesOf(columns);
    const TextEncoding encoding = buildTextEncoding(columns, columnNames);

    std::vector<veil::TypeFull> columnTypes;
    columnTypes.reserve(static_cast<size_t>(columns.size()));
    for (R_xlen_t i = 0; i < columns.size(); ++i)
    {
        const std::string name = static_cast<std::string>(columnNames[i]);
        columnTypes.push_back(rColumnType(VECTOR_ELT(SEXP(columns), i), name, encoding.levels()));
    }

    veil::Tree tree;
    veil::ObjStore objs;
    const veil::NodeId root = ingest(tree, objs, node, columnNames);
    tree.setRoot(root);
    veil::passAnnotateTypes(tree, columnTypes, objs);
    veil::passEncodeText(tree, encoding.mapping);
    veil::passInsertCoercions(tree);
    const std::vector<char> timeVarying = veil::passTagTimeVarying(tree);
    veil::passCheckTimeInvariance(tree, timeVarying);

    return dumpTree(tree, timeVarying, encoding.mapping);
}

// Reads a named list of columns into a ColumnSet and reports the constraints the scan reads off each
// -- NAs, whether the values are constant, and the per-type range -- so the column scan can be tested
// from R before the folding pass that consumes it exists.
[[cpp11::register]]
cpp11::list cpp_veil_scan(cpp11::list columns)
{
    const cpp11::strings columnNames = columnNamesOf(columns);
    const TextEncoding encoding = buildTextEncoding(columns, columnNames);

    veil::ColumnSet set;
    for (R_xlen_t i = 0; i < columns.size(); ++i)
    {
        const std::string name = static_cast<std::string>(columnNames[i]);
        SEXP column = VECTOR_ELT(SEXP(columns), i);
        readColumn(set, name, column, rColumnType(column, name, encoding.levels()),
                   encoding.indicesFor(static_cast<size_t>(i)));
    }

    cpp11::writable::strings names;
    cpp11::writable::strings types;
    cpp11::writable::logicals hasNAs;
    cpp11::writable::logicals constant;
    cpp11::writable::doubles mins;
    cpp11::writable::doubles maxes;
    cpp11::writable::logicals allIntegral;

    for (veil::ColumnId id = 0; id < static_cast<veil::ColumnId>(set.size()); ++id)
    {
        const veil::TypeWithConstraints scanned = veil::scanColumn(set.at(id));

        double low = NA_REAL;
        double high = NA_REAL;
        bool integral = true;
        if (const auto* d = std::get_if<veil::DoubleConstraint>(&scanned.specific))
        {
            if (scanned.hasValues) { low = d->min; high = d->max; }
            integral = d->allIntegral;
        }
        else if (const auto* n = std::get_if<veil::IntConstraint>(&scanned.specific))
        {
            if (scanned.hasValues) { low = static_cast<double>(n->min); high = static_cast<double>(n->max); }
        }
        else if (const auto* b = std::get_if<veil::BoolConstraint>(&scanned.specific))
        {
            if (scanned.hasValues) { low = b->min ? 1.0 : 0.0; high = b->max ? 1.0 : 0.0; }
        }

        names.push_back(set.name(id));
        types.push_back(typeName(scanned.type));
        hasNAs.push_back(static_cast<cpp11::r_bool>(scanned.hasNAs));
        constant.push_back(static_cast<cpp11::r_bool>(scanned.valuesAreConstant));
        mins.push_back(low);
        maxes.push_back(high);
        allIntegral.push_back(static_cast<cpp11::r_bool>(integral));
    }

    return cpp11::writable::list({
        "names"_nm = names,
        "types"_nm = types,
        "has_nas"_nm = hasNAs,
        "constant"_nm = constant,
        "min"_nm = mins,
        "max"_nm = maxes,
        "all_integral"_nm = allIntegral,
    });
}

// What the tree pipeline leaves behind for a caller to report.
struct PipelineResult final
{
    std::vector<char> timeVarying;
    std::vector<veil::Interval> intervals;
    size_t shared = 0; // How many nodes sharing merged away, so a test can see that it fired.
};

// Runs every tree pass, in order, over a tree that has already been ingested and given a root. One
// definition shared by the entry points below, so the pipeline cannot drift between them.
PipelineResult runTreePipeline(
    veil::Tree& tree,
    veil::ObjStore& objs,
    const std::vector<veil::TypeFull>& columnTypes,
    const std::vector<std::optional<veil::TypeWithConstraints>>& constraints,
    const veil::StringMapping& mapping,
    // What `.t` can be, in clicks. Defaulted to the representable calendar for the diagnostic entry
    // points, which are handed an expression and columns and no exposure to bound it with.
    veil::Interval timeInterval = veil::Interval::bounds(
        static_cast<double>(ValidDateStartClicks),
        static_cast<double>(ValidDateEndClicks)))
{
    veil::passAnnotateTypes(tree, columnTypes, objs);

    // Straight after annotation, which is what makes the operand types available, and before anything
    // that reasons about values: from here on no node is text, so the folding and interval passes
    // below never have to know that a category was ever written as a string.
    veil::passEncodeText(tree, mapping);

    veil::passInsertCoercions(tree);

    const std::vector<char> tagged = veil::passTagTimeVarying(tree);
    veil::passCheckTimeInvariance(tree, tagged);

    veil::passFoldConstants(tree, constraints);
    veil::passNarrowComparisons(tree);

    std::vector<veil::Interval> intervals =
        veil::passPropagateIntervals(tree, constraints, objs, timeInterval);
    veil::passFoldIntervalComparisons(tree, intervals);

    // Hoisting needs tags that cover every node, and narrowing has appended some since they were
    // taken. Tagging is non-mutating and cheap, so it is simply re-run rather than maintained.
    veil::passHoistFromIntegrate(tree, veil::passTagTimeVarying(tree));

    // Sharing runs LAST of the rewriting passes, because every one above it may edit a node in place
    // and doing that to a node with two parents corrupts the other. From here the tree is a DAG, and
    // anything added after this must build a new node rather than change an existing one.
    size_t shared = veil::passShareCommonSubtrees(tree);

    // `x * x` for an indicator x becomes x, which needs sharing to have run first -- a logical weight
    // reaches arithmetic as two separate to_double nodes, one per argument, and only sharing makes
    // them the same node so that a square is recognisable as one. Folding it then leaves an AEV's V
    // identical to its E, which is a redundancy that did not exist when sharing last ran, so sharing
    // runs again. That is the concrete reason CSE is worth running more than once.
    if (veil::passFoldIndicatorSquares(tree, constraints) > 0)
    {
        shared += veil::passShareCommonSubtrees(tree);
    }

    // The passes above append nodes -- narrowing adds a click threshold literal -- so the tags taken
    // before them no longer cover every node. Tagging is non-mutating and safe to repeat, so it is
    // re-run here against the tree as it finally stands.
    std::vector<char> timeVarying = veil::passTagTimeVarying(tree);

    // The intervals handed back are the ones the fold consumed, so they describe each node as it was
    // BEFORE folding. A comparison that folded to a literal still reports the 0-to-1 interval it had
    // as a comparison. That is only a reporting wrinkle on a debugging entry point, but it is why
    // these are not re-derived here.
    return PipelineResult{std::move(timeVarying), std::move(intervals), shared};
}

// What the block pipeline leaves behind for a caller to report.
struct BlockPipelineResult final
{
    size_t vectorOperands = 0; // How many time vectors the block names.
    size_t buffers = 0;        // How many it needs at once, which is what allocation is worth.
    size_t deathOnly = 0;      // How many are read at the death slot alone, and nowhere else.
    size_t slotEvaluations = 0; // Slots actually filled, over every instruction and individual.
    size_t dead = 0;            // Instructions writing a result nothing reads. Lowering emits none.
};

// Runs the TAC passes over a freshly lowered block. The tree phase decides WHAT to compute and has
// finished by the time this runs; these passes decide HOW it is laid out, and none of them changes
// an answer.
//
// One definition shared by every entry point, for the same reason runTreePipeline is: a block that
// went through a different set of passes depending on which entry point built it would be a
// difference nobody could see from the outside until an answer disagreed.
BlockPipelineResult runBlockPipeline(veil::Block& block)
{
    const std::vector<veil::VectorLiveRange> ranges = veil::passComputeLiveness(block);

    // Where each time vector is read comes first, because it says how much of the grid an
    // instruction has to fill; where it is stored comes after, and does not depend on it.
    const std::vector<veil::SlotDemand> demand = veil::passComputeSlotDemand(block);

    BlockPipelineResult layout;
    for (size_t id = 0; id < ranges.size(); ++id)
    {
        if (!ranges[id].defined) { continue; }
        ++layout.vectorOperands;
        if (demand[id].death && !demand[id].midpoints) { ++layout.deathOnly; }
    }

    layout.dead = veil::passFindDeadInstructions(block).count();

    block.assignSlotDemand(demand);
    layout.buffers = veil::passAllocateVectorBuffers(block, ranges);
    return layout;
}

// The record-chunking rule, exposed so it can be swept from R.
//
// The rule decides the order a total is accumulated in, so it is part of what an answer means -- and
// the properties that make it sound (every record in exactly one chunk, in order, with sizes differing
// by at most one) are worth checking over many record counts rather than at the handful of points the
// static assertions in RecordChunk.hpp can cover.
[[cpp11::register]]
cpp11::list cpp_veil_record_chunks(int records)
{
    if (records < 0) { cpp11::stop("A record count cannot be negative."); }

    const size_t count = static_cast<size_t>(records);
    cpp11::writable::integers startIndex;
    cpp11::writable::integers endIndex;
    for (size_t chunk = 0; chunk < veil::chunkCount(count); ++chunk)
    {
        const veil::RecordChunk range = veil::chunkOf(count, chunk);
        startIndex.push_back(static_cast<int>(range.startIndex));
        endIndex.push_back(static_cast<int>(range.endIndex));
    }

    return cpp11::writable::list({
        "start_index"_nm = startIndex,
        "end_index"_nm = endIndex,
        "records_per_chunk"_nm = static_cast<int>(veil::RecordsPerChunk),
    });
}

// The thread count a front end should ask for when the caller has not said.
//
// EXPOSED FOR THE SAME REASON AS THE TIME SCALES: the cap is deliberately low because CRAN checks on
// two cores, and it lives in ThreadPool.hpp so that every front end inherits it. A copy of the number
// in R would be a copy that stops matching.
[[cpp11::register]]
int cpp_veil_default_threads()
{
    return static_cast<int>(veil::DefaultThreadCount);
}

// The permitted time scales, in clicks.
//
// EXPOSED SO THAT R VALIDATES AGAINST THE ENGINE'S OWN LIST rather than a copy of it. A second
// spelling of this set in R would be a place for the two to drift, and the symptom would be either a
// scale R accepts and the engine refuses, or -- worse -- one R refuses that would have been fine.
[[cpp11::register]]
cpp11::integers cpp_veil_time_scales()
{
    cpp11::writable::integers scales;
    for (const int clicks : veil::timeScaleClickOptions) { scales.push_back(clicks); }
    return scales;
}

// The canonical clicks-to-years conversion, exposed so its exactness can be swept from R. What a
// click and a year mean to each other is a property the whole system depends on agreeing about, and
// the exactness it has at whole years, ages, months and quarters follows from ClicksPerYear's
// particular value rather than from the arithmetic, so it is worth a test that fires if that value
// ever changes.
[[cpp11::register]]
cpp11::doubles cpp_veil_to_years(cpp11::integers clicks)
{
    cpp11::writable::doubles years;
    for (R_xlen_t i = 0; i < clicks.size(); ++i)
    {
        years.push_back(clicks[i] == NA_INTEGER ? NA_REAL : yearsFromClicks(clicks[i]));
    }
    return years;
}

// The reverse conversion, per the datey specification: scale, round half to even, then take the
// integer. A value whose click count falls outside `int` reports NA rather than overflowing.
[[cpp11::register]]
cpp11::integers cpp_veil_to_clicks(cpp11::doubles years)
{
    cpp11::writable::integers clicks;
    for (R_xlen_t i = 0; i < years.size(); ++i)
    {
        const double value = years[i];
        const double scaled = ISNA(value) || !R_finite(value) ? NA_REAL : clicksFromYears(value);
        const bool representable = !ISNA(scaled)
            && scaled >= static_cast<double>(std::numeric_limits<int>::min())
            && scaled <= static_cast<double>(std::numeric_limits<int>::max());
        clicks.push_back(representable ? static_cast<int>(scaled) : NA_INTEGER);
    }
    return clicks;
}

// Runs the whole tree pipeline and dumps the rewritten tree, so the folding can be tested from R.
[[cpp11::register]]
cpp11::list cpp_veil_fold(cpp11::list node, cpp11::list columns)
{
    const cpp11::strings columnNames = columnNamesOf(columns);
    std::vector<veil::TypeFull> columnTypes;
    std::vector<std::optional<veil::TypeWithConstraints>> constraints;
    veil::ColumnSet set;
    const TextEncoding encoding = prepareColumns(columns, columnNames, set, columnTypes, constraints);

    veil::Tree tree;
    veil::ObjStore objs;
    tree.setRoot(ingest(tree, objs, node, columnNames));
    const PipelineResult result =
        runTreePipeline(tree, objs, columnTypes, constraints, encoding.mapping);

    return dumpTree(tree, result.timeVarying, encoding.mapping);
}

// The exposure, if the data carries one. This asks rather than insists: an entry point that does not
// need a time vector still benefits from supplying one when it is there, because an expression that
// turns out to sample time then fails on what is actually wrong with it rather than on a column the
// caller had no reason to think was needed.
std::optional<veil::ExposureColumns> findExposure(
    const cpp11::strings& columnNames,
    const std::vector<veil::TypeFull>& columnTypes)
{
    veil::ExposureColumns exposure;

    // Tracked one flag per column rather than by counting. R lists may carry duplicate names, and a
    // count of three is also what two `E2R_start` columns and one `E2R_end` would give -- which
    // would then describe an exposure with no end.
    bool haveStart = false;
    bool haveEnd = false;
    bool haveDied = false;

    for (R_xlen_t i = 0; i < columnNames.size(); ++i)
    {
        const std::string name = static_cast<std::string>(columnNames[i]);
        const veil::ColumnId id = static_cast<veil::ColumnId>(i);
        if (name == "E2R_start" && columnTypes[id].type == veil::Type::Datey) { exposure.start = id; haveStart = true; }
        else if (name == "E2R_end" && columnTypes[id].type == veil::Type::Datey) { exposure.end = id; haveEnd = true; }
        else if (name == "E2R_died" && columnTypes[id].type == veil::Type::Bool) { exposure.died = id; haveDied = true; }
    }
    if (!haveStart || !haveEnd || !haveDied) { return std::nullopt; }

    exposure.deltaTClicks = veil::defaultTimeScaleClicks;
    return exposure;
}

// Lines up the columns the reader took with the ColumnIds a field reference carries.
//
// They are two different numberings: a field's ColumnId indexes the caller's whole list, while the
// ColumnSet numbers only the columns it could read, in the order it read them. A text column in the
// middle of the list would otherwise shift everything after it. The entry is null where a column was
// not read, and a block that binds one of those is an error the interpreter raises by name.
std::vector<const veil::ColumnView*> viewsByColumnId(
    const veil::ColumnSet& set,
    const cpp11::strings& columnNames)
{
    std::vector<const veil::ColumnView*> views(static_cast<size_t>(columnNames.size()), nullptr);
    for (R_xlen_t i = 0; i < columnNames.size(); ++i)
    {
        if (const std::optional<veil::ColumnId> id = set.find(static_cast<std::string>(columnNames[i])))
        {
            views[static_cast<size_t>(i)] = &set.at(*id);
        }
    }
    return views;
}

// Runs a lowered block once per individual and hands back one R value per individual, in the R type
// the block's output operand calls for.
// Loads one individual for an entry point that wants a value for every one of them.
//
// These entry points accept no include, so nothing can empty an individual and `loadRecord` always
// answers true. Checked rather than assumed: were an include ever plumbed through to one of them, the
// alternative to this complaint is a register holding the PREVIOUS individual's answer.
void loadEveryRecord(
    veil::Interpreter& interpreter,
    const std::vector<const veil::ColumnView*>& views,
    size_t record)
{
    if (!interpreter.loadRecord(views, record))
    {
        cpp11::stop("veil: an individual was excluded by an entry point that has no include.");
    }
}

cpp11::sexp runOverRecords(
    const veil::Block& block,
    const std::vector<const veil::ColumnView*>& views,
    R_xlen_t records)
{
    veil::Interpreter interpreter(block);
    const veil::TypeFull outputType = block.operandAt(block.outputs().front()).type;

    switch (outputType.type)
    {
        case veil::Type::Double:
        {
            cpp11::writable::doubles out(records);
            for (R_xlen_t i = 0; i < records; ++i)
            {
                loadEveryRecord(interpreter, views, static_cast<size_t>(i));
                interpreter.run();
                out[i] = interpreter.outputAsNumber(0);
            }
            return cpp11::sexp(SEXP(out));
        }
        case veil::Type::Bool:
        {
            cpp11::writable::logicals out(records);
            for (R_xlen_t i = 0; i < records; ++i)
            {
                loadEveryRecord(interpreter, views, static_cast<size_t>(i));
                interpreter.run();
                out[i] = static_cast<cpp11::r_bool>(interpreter.outputAsLogical(0));
            }
            return cpp11::sexp(SEXP(out));
        }
        case veil::Type::Datey:
        case veil::Type::Durationy:
        {
            // Bare clicks, with no datey class attached: this is a test and diagnostic entry point,
            // and a test comparing against R does so through unclass().
            cpp11::writable::integers out(records);
            for (R_xlen_t i = 0; i < records; ++i)
            {
                loadEveryRecord(interpreter, views, static_cast<size_t>(i));
                interpreter.run();
                out[i] = interpreter.outputAsClicks(0);
            }
            return cpp11::sexp(SEXP(out));
        }
        default:
            break;
    }

    cpp11::stop("veil cannot yet hand back a result of type '%s'.", typeName(outputType));
}

// Compiles SEVERAL expressions as one calculation and evaluates them all together, which is the shape
// every recipe has: an AEV wants three values out of one pass over the data, a log-likelihood's second
// differential wants one per pair of parameters.
//
// The roots share one arena, so a sub-expression two of them have in common is lowered ONCE. That is
// what the instruction count reports, and it is the property worth testing -- the answers alone would
// look the same whether or not anything was shared.
//
// Values come back as doubles whatever each output's type, because the point here is the sharing
// rather than the marshalling; a datey output arrives as its clicks. `types` says what each one is.
[[cpp11::register]]
cpp11::list cpp_veil_eval_multi(cpp11::list nodes, cpp11::list columns)
{
    const cpp11::strings columnNames = columnNamesOf(columns);
    std::vector<veil::TypeFull> columnTypes;
    std::vector<std::optional<veil::TypeWithConstraints>> constraints;
    veil::ColumnSet set;
    const TextEncoding encoding = prepareColumns(columns, columnNames, set, columnTypes, constraints);

    veil::Tree tree;
    veil::ObjStore objs;
    for (R_xlen_t i = 0; i < nodes.size(); ++i)
    {
        tree.addRoot(ingest(tree, objs, cpp11::as_cpp<cpp11::list>(nodes[i]), columnNames));
    }
    if (tree.roots().empty()) { cpp11::stop("veil needs at least one expression to evaluate."); }

    const PipelineResult result =
        runTreePipeline(tree, objs, columnTypes, constraints, encoding.mapping);
    veil::Block block =
        veil::passLowerToBlock(
            tree, objs, result.timeVarying,
            findExposure(columnNames, columnTypes));
    const BlockPipelineResult layout = runBlockPipeline(block);

    const R_xlen_t records = columns.size() == 0 ? 0 : Rf_xlength(VECTOR_ELT(SEXP(columns), 0));
    const std::vector<const veil::ColumnView*> views = viewsByColumnId(set, columnNames);
    const size_t outputCount = block.outputs().size();

    std::vector<cpp11::writable::doubles> collected(outputCount, cpp11::writable::doubles(records));
    veil::Interpreter interpreter(block);
    for (R_xlen_t i = 0; i < records; ++i)
    {
        loadEveryRecord(interpreter, views, static_cast<size_t>(i));
        interpreter.run();
        for (size_t j = 0; j < outputCount; ++j)
        {
            const veil::TypeFull type = block.operandAt(block.outputs()[j]).type;
            switch (type.type)
            {
                case veil::Type::Double: collected[j][i] = interpreter.outputAsNumber(j); break;
                case veil::Type::Bool:
                    collected[j][i] = interpreter.outputAsLogical(j) ? 1.0 : 0.0;
                    break;
                case veil::Type::Datey:
                case veil::Type::Durationy:
                    collected[j][i] = static_cast<double>(interpreter.outputAsClicks(j));
                    break;
                default: cpp11::stop("veil cannot yet hand back a result of type '%s'.", typeName(type));
            }
        }
    }

    cpp11::writable::list values;
    cpp11::writable::strings types;
    for (size_t j = 0; j < outputCount; ++j)
    {
        values.push_back(collected[j]);
        types.push_back(typeName(block.operandAt(block.outputs()[j]).type));
    }

    cpp11::writable::strings monikers;
    for (const veil::Instruction& instruction : block.body())
    {
        monikers.push_back(std::string(veil::moniker(instruction.op)));
    }

    return cpp11::writable::list({
        "values"_nm = values,
        "types"_nm = types,
        "monikers"_nm = monikers,
        "output_count"_nm = static_cast<int>(outputCount),
        "instruction_count"_nm = static_cast<int>(block.body().size()),
        "operand_count"_nm = static_cast<int>(block.operandCount()),
        "vector_operand_count"_nm = static_cast<int>(layout.vectorOperands),
        "buffer_count"_nm = static_cast<int>(layout.buffers),
        "death_only_count"_nm = static_cast<int>(layout.deathOnly),
        "dead_instruction_count"_nm = static_cast<int>(layout.dead),
        "shared_nodes"_nm = static_cast<int>(result.shared),
        "node_count"_nm = static_cast<int>(tree.size()),
    });
}

// Compiles an expression the whole way and evaluates it for every individual, so that what veil
// computes can be checked against what R computes for the same expression over the same columns.
//
// It covers the scalar spine: a time-invariant value per individual. Anything needing the time
// vector -- `.t`, a mortality, an integration -- is refused by lowering with a message saying so.
//
// The report alongside the values is the lowered block, which is what makes the compilation itself
// testable: an expression the tree phase settled to a constant lowers to no instructions at all, and
// a field read twice lowers to one column-bound operand rather than two.
//
// Note that every individual gets a value, including for an expression that mentions no field at
// all: the block runs once per row, where R's own evaluation of a constant expression would give a
// single value.
[[cpp11::register]]
cpp11::list cpp_veil_eval(cpp11::list node, cpp11::list columns)
{
    const cpp11::strings columnNames = columnNamesOf(columns);
    std::vector<veil::TypeFull> columnTypes;
    std::vector<std::optional<veil::TypeWithConstraints>> constraints;
    veil::ColumnSet set;
    const TextEncoding encoding = prepareColumns(columns, columnNames, set, columnTypes, constraints);

    veil::Tree tree;
    veil::ObjStore objs;
    tree.setRoot(ingest(tree, objs, node, columnNames));
    const PipelineResult result =
        runTreePipeline(tree, objs, columnTypes, constraints, encoding.mapping);

    // The exposure is handed over when the data has one, even though this entry point wants a single
    // value per individual. An expression that samples time is then refused for having no such
    // value, which is the useful complaint, rather than for a column the caller never mentioned.
    veil::Block block =
        veil::passLowerToBlock(
            tree, objs, result.timeVarying,
            findExposure(columnNames, columnTypes));
    if (block.outputs().size() != 1)
    {
        cpp11::stop("A veil expression lowers to exactly one output.");
    }
    runBlockPipeline(block);

    // The individual count comes from the data as R passed it, not from the ColumnSet: an expression
    // that reads no column still has one value per individual.
    const R_xlen_t records = columns.size() == 0
        ? 0 : Rf_xlength(VECTOR_ELT(SEXP(columns), 0));

    const cpp11::sexp values = runOverRecords(block, viewsByColumnId(set, columnNames), records);

    cpp11::writable::strings monikers;
    for (const veil::Instruction& instruction : block.body())
    {
        monikers.push_back(std::string(veil::moniker(instruction.op)));
    }

    return cpp11::writable::list({
        "values"_nm = values,
        "type"_nm = typeName(block.operandAt(block.outputs().front()).type),
        "monikers"_nm = monikers,
        "operand_count"_nm = static_cast<int>(block.operandCount()),
        "constant_count"_nm = static_cast<int>(block.constants().size()),
        "column_count"_nm = static_cast<int>(block.columns().size()),
    });
}

// Resolves one of the exposure columns, and checks it is the type the time grid needs. A clearer
// failure than letting the field resolver complain about a pronoun the user never wrote.
veil::ColumnId exposureColumn(
    const cpp11::strings& columnNames,
    const std::vector<veil::TypeFull>& columnTypes,
    const char* name,
    veil::Type wanted)
{
    for (R_xlen_t i = 0; i < columnNames.size(); ++i)
    {
        if (static_cast<std::string>(columnNames[i]) == name)
        {
            const veil::ColumnId id = static_cast<veil::ColumnId>(i);
            if (columnTypes[id].type != wanted)
            {
                cpp11::stop(
                    "Column `%s` is a %s, which the time vector cannot use.",
                    name,
                    typeName(columnTypes[id]));
            }
            return id;
        }
    }
    cpp11::stop("Integrating needs an `%s` column, which the data does not have.", name);
}

// Compiles `node` wrapped in a finalising op and runs it for every individual.
//
// The wrapper goes on BEFORE the tree pipeline rather than after, so that every pass sees the whole
// expression. Wrapping afterwards would leave the tags and the intervals describing a tree that no
// longer exists, which is the failure mode the time-varying tags already had once.
cpp11::writable::doubles runFinalised(
    cpp11::list node,
    const cpp11::strings& columnNames,
    const std::vector<veil::TypeFull>& columnTypes,
    const std::vector<std::optional<veil::TypeWithConstraints>>& constraints,
    const veil::StringMapping& mapping,
    const std::vector<const veil::ColumnView*>& views,
    R_xlen_t records,
    veil::Op finaliseOp,
    const veil::ExposureColumns& exposure,
    SEXP include,
    cpp11::writable::integers* slotCounts,
    cpp11::writable::strings* monikers,
    BlockPipelineResult* layout = nullptr)
{
    veil::Tree tree;
    veil::ObjStore objs;
    const veil::NodeId integrand = ingest(tree, objs, node, columnNames);
    tree.setRoot(tree.buildCall(finaliseOp, {integrand}));

    // The include is read BEFORE the pipeline runs, because its gates are ingested into this same
    // tree and have to be typed and coerced along with everything else. They are not reachable from
    // the root, which is fine: the passes that matter for a gate sweep every node rather than walking
    // down from the root. Interval propagation does walk from the root, so a gate simply carries no
    // interval and nothing folds inside it -- sound, if not yet clever.
    std::optional<veil::ObjId> includeObj;
    if (include != R_NilValue)
    {
        if (!Rf_inherits(include, "include"))
        {
            cpp11::stop("`include` must be an `include` object.");
        }
        includeObj = readInclude(tree, objs, include, columnNames);
    }

    const PipelineResult result = runTreePipeline(tree, objs, columnTypes, constraints, mapping);
    veil::Block block =
        veil::passLowerToBlock(tree, objs, result.timeVarying, exposure, includeObj);

    const BlockPipelineResult blockLayout = runBlockPipeline(block);
    if (layout != nullptr) { *layout = blockLayout; }

    if (monikers != nullptr)
    {
        for (const veil::Instruction& instruction : block.body())
        {
            monikers->push_back(std::string(veil::moniker(instruction.op)));
        }
    }

    veil::Interpreter interpreter(block);
    cpp11::writable::doubles out(records);
    for (R_xlen_t i = 0; i < records; ++i)
    {
        // This entry point DOES take an include, so an individual can be emptied by it. Their answer
        // is zero -- the integral of nothing, or the value at a death that falls outside the window
        // -- and it has to be written rather than left, because the register still holds whoever came
        // before them. Nothing else about them is computed.
        if (!interpreter.loadRecord(views, static_cast<size_t>(i)))
        {
            out[i] = 0.0;
            if (slotCounts != nullptr) { slotCounts->push_back(interpreter.slotCount()); }
            continue;
        }

        interpreter.run();
        out[i] = interpreter.outputAsNumber(0);
        if (slotCounts != nullptr) { slotCounts->push_back(interpreter.slotCount()); }
    }
    if (layout != nullptr) { layout->slotEvaluations = interpreter.slotEvaluations(); }
    return out;
}

// Integrates an expression over each individual's exposure, and takes its value at death, so that
// the time vector can be checked from R.
//
// The integral is a rectangular midpoint rule over intervals of one part in `intervals_per_year`,
// with a final interval of that length or less. It is EXACT for an integrand that is constant or
// linear in time, which is what makes it checkable against an answer worked out on paper rather than
// against a second implementation of the same rule.
//
// `integrate` and `died_value` are deliberately absent from the expression whitelist -- they are
// internal, and a user is not meant to depend on them -- so this entry point applies them rather
// than accepting them written out.
// `include` may be R's NULL, meaning the whole exposure counts. When it is an `include` object the
// exposure is clipped to its interval and the death flag is recomputed, per the specification:
// nu' = max(nu, from), tau' = min(tau, to), and an individual counts as a death only if the clip did
// not cut the end off. An individual the include empties integrates to zero rather than erroring.
[[cpp11::register]]
cpp11::list cpp_veil_integrate(
    cpp11::list node,
    cpp11::list columns,
    int time_scale,
    SEXP include)
{
    const cpp11::strings columnNames = columnNamesOf(columns);
    std::vector<veil::TypeFull> columnTypes;
    std::vector<std::optional<veil::TypeWithConstraints>> constraints;
    veil::ColumnSet set;
    const TextEncoding encoding = prepareColumns(columns, columnNames, set, columnTypes, constraints);

    veil::ExposureColumns exposure;
    exposure.start = exposureColumn(columnNames, columnTypes, "E2R_start", veil::Type::Datey);
    exposure.end = exposureColumn(columnNames, columnTypes, "E2R_end", veil::Type::Datey);
    exposure.died = exposureColumn(columnNames, columnTypes, "E2R_died", veil::Type::Bool);
    exposure.deltaTClicks = veil::validateTimeScaleClicks(time_scale);

    const std::vector<const veil::ColumnView*> views = viewsByColumnId(set, columnNames);
    const R_xlen_t records = columns.size() == 0 ? 0 : Rf_xlength(VECTOR_ELT(SEXP(columns), 0));

    cpp11::writable::integers slotCounts;
    cpp11::writable::strings monikers;
    BlockPipelineResult layout;
    const cpp11::writable::doubles integral =
        runFinalised(
            node, columnNames, columnTypes, constraints, encoding.mapping, views, records,
            veil::Op::Integrate, exposure, include, &slotCounts, &monikers, &layout);
    const cpp11::writable::doubles atDeath =
        runFinalised(
            node, columnNames, columnTypes, constraints, encoding.mapping, views, records,
            veil::Op::DiedValue, exposure, include, nullptr, nullptr);

    // The layout reported is the integral's, since that is the block with the time vector in it.
    return cpp11::writable::list({
        "integral"_nm = integral,
        "died_value"_nm = atDeath,
        "slot_count"_nm = slotCounts,
        "monikers"_nm = monikers,
        "delta_t_clicks"_nm = exposure.deltaTClicks,
        "vector_operand_count"_nm = static_cast<int>(layout.vectorOperands),
        "buffer_count"_nm = static_cast<int>(layout.buffers),
        "death_only_count"_nm = static_cast<int>(layout.deathOnly),
        "dead_instruction_count"_nm = static_cast<int>(layout.dead),
        "slot_evaluations"_nm = static_cast<int>(layout.slotEvaluations),
    });
}

// The A/E/V calculation: one pass over the data producing one A, one E and one V.
//
// `mortality` and `weight` are pronoun ASTs -- a mortality object crosses as an obj leaf like any
// other -- and `include` may be R's NULL. The weight may be NULL too, meaning a weight of one, which
// makes A/E/V a count of lives.
//
// `A`, `E` and `V` ARE THE TOTALS, added up in the core. The summation is part of the calculation
// rather than something the host does afterwards: it has to mean the same thing through every front
// end, and its ORDER has to be a function of the data alone or a threaded run would disagree with
// itself. See Engine.hpp and RecordChunk.hpp.
//
// `contributions` holds the per-individual values as well, and is DIAGNOSTIC. It is what lets the
// analytic oracles be written per individual, and it is the oracle for the accumulation itself, since
// a plain sum of it must agree with the chunked total. It goes when the batch entry point lands.
namespace
{

// Runs R's interrupt check inside a top-level context so that its longjmp can be caught.
void checkInterruptImmediately(void*) { R_CheckUserInterrupt(); }

// Answers whether the user has asked R to stop, WITHOUT longjmping.
//
// `R_CheckUserInterrupt` jumps to the top level when an interrupt is pending, which from inside a
// calculation would skip every destructor between here and there -- and from a worker thread would
// take the process down. `R_ToplevelExec` runs a function in a fresh top-level context and answers
// FALSE if it jumped, which is how the jump becomes a bool.
//
// IT CONSUMES THE INTERRUPT: the condition is caught by that context rather than delivered, so the
// pending flag is cleared and the caller must raise the error itself. That is exactly the order the
// engine needs -- stop the workers, join them, and only then raise -- rather than a constraint to
// work around.
//
// MAIN THREAD ONLY, which runInParallel guarantees: it calls this from the thread that called it and
// from no other. Nothing about R crosses into the core; the core takes a `bool()` and this is one.
bool userInterruptIsPending()
{
    return R_ToplevelExec(checkInterruptImmediately, nullptr) == FALSE;
}

} // namespace

// What building one AEV specification leaves behind: the block to run, and what the pipeline did on
// the way, which the entry points report so a test can see that a pass fired.
struct AevBuild final
{
    veil::Block block;
    size_t shared = 0;
    BlockPipelineResult layout;
};

// What `.t` can be over a whole dataset, in clicks: the earliest exposure start to the latest
// exposure end.
//
// SOUND BECAUSE OF WHERE THE INTEGRATION SAMPLES. Every sample point of an individual's time vector
// lies within their own exposure, and the death slot is its end exactly; an include only narrows
// that further. So no `.t` anywhere in the crossing falls outside these two columns' extremes.
//
// A scan that gives nothing away -- a column with missing values, which the reader refuses anyway,
// or one with no rows -- leaves the bound where it was, so this can only tighten and never widens.
veil::Interval exposureTimeInterval(
    const veil::ExposureColumns& exposure,
    const std::vector<std::optional<veil::TypeWithConstraints>>& constraints)
{
    const veil::Interval whole = veil::Interval::bounds(
        static_cast<double>(ValidDateStartClicks),
        static_cast<double>(ValidDateEndClicks));

    const auto scanned = [&constraints](veil::ColumnId column)
    {
        return column < constraints.size() && constraints[column].has_value()
            ? veil::detail::columnInterval(*constraints[column])
            : veil::Interval::unknown();
    };

    const veil::Interval start = scanned(exposure.start);
    const veil::Interval end = scanned(exposure.end);

    return veil::Interval::bounds(std::isfinite(start.lo) ? start.lo : whole.lo,
                                  std::isfinite(end.hi) ? end.hi : whole.hi);
}

// Compiles ONE AEV specification -- a log-mu expression, an optional weight and an optional include --
// into a block ready to run.
//
// ONE DEFINITION SHARED BY THE ENTRY POINTS, for the reason given on runTreePipeline: a block built a
// slightly different way depending on which entry point asked for it is a difference nobody can see
// from the outside until two answers disagree. `cpp_veil_aev` compiles one of these and
// `cpp_veil_run` compiles several, and that is the only difference between them.
//
// The dataset-wide arguments are passed in already prepared, because every spec in a batch shares
// them: reading the columns once is the point of a batch entry point.
AevBuild buildAev(
    cpp11::list mortality,
    SEXP weight,
    SEXP valSimilarity,
    SEXP valDistance,
    SEXP include,
    const cpp11::strings& columnNames,
    const std::vector<veil::TypeFull>& columnTypes,
    const std::vector<std::optional<veil::TypeWithConstraints>>& constraints,
    const veil::StringMapping& mapping,
    const veil::ExposureColumns& exposure)
{
    veil::Tree tree;
    veil::ObjStore objs;

    const veil::NodeId logMu = ingest(tree, objs, mortality, columnNames);
    const veil::NodeId weightNode = weight == R_NilValue
        ? tree.buildLitDouble(1.0)
        : ingest(tree, objs, cpp11::as_cpp<cpp11::list>(weight), columnNames);

    // The two spellings of the second weighting factor are one quantity, so only one may be given.
    // Which one was written travels all the way to the recipe rather than being converted here; see
    // `SimilarityForm` for why.
    if (valSimilarity != R_NilValue && valDistance != R_NilValue)
    {
        cpp11::stop("Give either `val_similarity` or `val_distance`, not both.");
    }

    const SEXP secondFactor = valSimilarity != R_NilValue ? valSimilarity : valDistance;
    const veil::SimilarityForm form = valSimilarity != R_NilValue
        ? veil::SimilarityForm::Similarity
        : veil::SimilarityForm::Distance;
    const veil::NodeId similarityNode = secondFactor == R_NilValue
        ? veil::invalidNodeId
        : ingest(tree, objs, cpp11::as_cpp<cpp11::list>(secondFactor), columnNames);

    const veil::AevRoots roots = veil::buildAevRecipe(tree, logMu, weightNode, similarityNode, form);
    tree.addRoot(roots.a);
    tree.addRoot(roots.e);
    tree.addRoot(roots.v);

    // Read before the pipeline, so an include's gates are typed and coerced with everything else.
    std::optional<veil::ObjId> includeObj;
    if (include != R_NilValue)
    {
        if (!Rf_inherits(include, "include")) { cpp11::stop("`include` must be an `include` object."); }
        includeObj = readInclude(tree, objs, include, columnNames);
    }

    const PipelineResult result = runTreePipeline(tree, objs, columnTypes, constraints, mapping,
                                                  exposureTimeInterval(exposure, constraints));

    // Said in the user's own words, since they wrote one spelling and not the other. Only a
    // violation that is certain before any record is read reaches here; the pass says why it
    // refuses nothing weaker.
    switch (veil::passCheckSimilarityRange(tree, result.intervals, similarityNode, form))
    {
        case veil::SimilarityViolation::Above:
            if (form == veil::SimilarityForm::Distance)
            {
                cpp11::stop("`val_distance` cannot be negative, and this one always is.");
            }
            cpp11::stop("`val_similarity` must lie between 0 and 1, and this one is always above 1.");
            break;
        case veil::SimilarityViolation::Below:
            cpp11::stop("`val_similarity` must lie between 0 and 1, and this one is always below 0.");
            break;
        case veil::SimilarityViolation::None:
            break;
    }

    veil::Block block =
        veil::passLowerToBlock(tree, objs, result.timeVarying, exposure, includeObj);
    const BlockPipelineResult layout = runBlockPipeline(block);

    if (block.outputs().size() != 3) { cpp11::stop("An AEV lowers to exactly three outputs."); }

    return AevBuild{std::move(block), result.shared, layout};
}

// The overdispersion a calculation was given, checked.
//
// MANDATORY AT EVERY ENTRY POINT, WITH NO DEFAULT ANYWHERE (Tim, 2026-08-06). B060 opens by saying
// that failing to allow for overdispersion underestimates uncertainty and overfits models, so a
// default would be exactly an invitation to ignore it. The comparison is written so that a NaN is
// refused along with zero and the negatives.
inline double checkedOverdispersion(double overdispersion)
{
    if (!(overdispersion > 0.0)) { cpp11::stop("`overdispersion` must be a positive number."); }
    return overdispersion;
}

// V IS A VARIANCE, NOT A SECOND MOMENT (Tim, 2026-07-31: "V should always be correct"). The engine
// integrates Ew^2; the overdispersion supplied for this calculation turns that into
// Var(Aw - Ew) = Omega E Ew^2, which is what an `aev`'s confidence intervals and residuals read. A V
// travelling without its own overdispersion can meet the wrong one downstream, which is an error
// rather than an inconvenience.
//
// APPLIED ON THE WAY OUT, ONCE, AND NOWHERE ELSE. Scaling the output slot rather than the integrand
// is what keeps passFoldIndicatorSquares and the shared-operand collapse firing, so w^2 = w still
// folds and V = E survives on the INTEGRALS -- it is only at the output that V = Omega E for an
// indicator weight. It also keeps the scaling clear of the chunked fold, so threading still cannot
// move a number.
inline double varianceFromSecondMoment(double secondMoment, double overdispersion)
{
    return secondMoment * overdispersion;
}

[[cpp11::register]]
cpp11::list cpp_veil_aev(
    cpp11::list mortality,
    SEXP weight,
    cpp11::list columns,
    int time_scale,
    SEXP include,
    double overdispersion,
    // Spelled at every call site because cpp11 does not carry a C++ default into the R wrapper it
    // generates. That is no loss here: it makes the calls that pin the single-threaded path say so.
    int threads)
{
    const cpp11::strings columnNames = columnNamesOf(columns);
    std::vector<veil::TypeFull> columnTypes;
    std::vector<std::optional<veil::TypeWithConstraints>> constraints;
    veil::ColumnSet set;
    const TextEncoding encoding = prepareColumns(columns, columnNames, set, columnTypes, constraints);

    veil::ExposureColumns exposure;
    exposure.start = exposureColumn(columnNames, columnTypes, "E2R_start", veil::Type::Datey);
    exposure.end = exposureColumn(columnNames, columnTypes, "E2R_end", veil::Type::Datey);
    exposure.died = exposureColumn(columnNames, columnTypes, "E2R_died", veil::Type::Bool);
    exposure.deltaTClicks = veil::validateTimeScaleClicks(time_scale);

    const double dispersion = checkedOverdispersion(overdispersion);

    const AevBuild built =
        buildAev(mortality, weight, R_NilValue, R_NilValue, include, columnNames, columnTypes,
                 constraints, encoding.mapping, exposure);
    const veil::Block& block = built.block;
    const BlockPipelineResult& layout = built.layout;

    const R_xlen_t records = columns.size() == 0 ? 0 : Rf_xlength(VECTOR_ELT(SEXP(columns), 0));
    const std::vector<const veil::ColumnView*> views = viewsByColumnId(set, columnNames);

    if (threads < 0) { cpp11::stop("`threads` cannot be negative."); }

    // THE TOTALS ARE THE ANSWER, and the engine computes them. The contributions come back alongside
    // for the tests only -- see CalculationResult -- and go when the batch entry point lands.
    //
    // THREADING CANNOT MOVE A NUMBER, which is the only reason a thread count is safe to expose at
    // all: the chunk partials are folded in chunk order whoever ran them. `threads` defaults to 1 so
    // that every caller that predates threading keeps exactly the path it had; the user-facing
    // default is veil::DefaultThreadCount, and it belongs on `aev()` when that lands.
    const size_t threadCount = static_cast<size_t>(threads);
    const veil::CalculationResult calculation = veil::runCalculation(
        block, views, static_cast<size_t>(records), true, threadCount, userInterruptIsPending);

    // Raised here rather than delivered by R, because the poll consumed the pending interrupt. Every
    // worker has been joined by the time this runs, so there is nothing still touching the results.
    if (calculation.interrupted) { cpp11::stop("The veil calculation was interrupted."); }

    cpp11::writable::doubles a(records);
    cpp11::writable::doubles e(records);
    cpp11::writable::doubles v(records);
    for (R_xlen_t i = 0; i < records; ++i)
    {
        const size_t base = static_cast<size_t>(i) * 3;
        a[i] = calculation.contributions[base];
        e[i] = calculation.contributions[base + 1];
        v[i] = varianceFromSecondMoment(calculation.contributions[base + 2], dispersion);
    }

    cpp11::writable::strings monikers;
    for (const veil::Instruction& instruction : block.body())
    {
        monikers.push_back(std::string(veil::moniker(instruction.op)));
    }

    return cpp11::writable::list({
        "A"_nm = calculation.totals[0],
        "E"_nm = calculation.totals[1],
        "V"_nm = varianceFromSecondMoment(calculation.totals[2], dispersion),
        "contributions"_nm = cpp11::writable::list({"A"_nm = a, "E"_nm = e, "V"_nm = v}),
        "chunk_count"_nm = static_cast<int>(calculation.chunks),
        "records_included"_nm = static_cast<int>(calculation.recordsIncluded),
        "monikers"_nm = monikers,
        "instruction_count"_nm = static_cast<int>(block.body().size()),
        "shared_nodes"_nm = static_cast<int>(built.shared),
        "vector_operand_count"_nm = static_cast<int>(layout.vectorOperands),
        "buffer_count"_nm = static_cast<int>(layout.buffers),
        "death_only_count"_nm = static_cast<int>(layout.deathOnly),
        "dead_instruction_count"_nm = static_cast<int>(layout.dead),
        "slot_evaluations"_nm = static_cast<int>(calculation.slotEvaluations),

        // THE ONLY OBSERVABLE THAT A THREAD COUNT WAS HONOURED. Nothing numeric can see threading --
        // that is the entire design -- so without this a test cannot tell four threads from a silent
        // fall back to one, and the bit-identity test would pass just as happily either way.
        "threads_used"_nm = static_cast<int>(veil::resolveThreadCount(threadCount)),
    });
}

namespace
{

// A named element of an R list, or NULL when the list has not got one. Looked up by hand rather than
// with cpp11's own name indexing, which raises for a missing name -- an absent `weight` or `include`
// is the ordinary case here, not an error.
SEXP elementOrNull(const cpp11::list& list, const char* name)
{
    const SEXP names = Rf_getAttrib(SEXP(list), R_NamesSymbol);
    if (names == R_NilValue) { return R_NilValue; }

    const R_xlen_t count = Rf_xlength(names);
    for (R_xlen_t i = 0; i < count; ++i)
    {
        if (std::string(CHAR(STRING_ELT(names, i))) == name) { return VECTOR_ELT(SEXP(list), i); }
    }
    return R_NilValue;
}

} // namespace

// Runs SEVERAL AEV specifications over ONE dataset in a single crossing.
//
// THIS IS WHERE OP-LEVEL PARALLELISM BEGINS. It is the primary axis -- a user runs several A/E/Vs and
// several fits together, so there are usually more operations than cores -- but until now every entry
// point compiled exactly one block, so there was never more than one operation to schedule and the
// record-chunk axis was carrying work it was only ever meant to supplement. `runCalculations`
// dispatches (block, chunk) pairs as one flat task list; this is what gives it several blocks.
//
// `specs` is a list of lists, each with a `mortality` and an `overdispersion`, and optionally a
// `weight`, an `include` and one of `val_similarity` / `val_distance`. The dataset, the time scale
// and the thread count are shared, which is the point: the columns are read and prepared once for
// the whole batch.
//
// OVERDISPERSION IS PER SPEC, the time scale is not. Overdispersion is a scalar on each block's V
// output slot, so entries in one batch may legitimately disagree about it -- which is what lets a
// user override it for one entry. The time scale is built into the one ExposureColumns every block
// shares, so it cannot vary within a crossing as things stand.
//
// EVERY SPEC IS COMPILED BEFORE ANY IS RUN. A batch with one bad specification then fails without
// having spent a dataset walk on the others, and it fails the same way whichever spec is bad.
//
// `keep_contributions` is diagnostic and OFF for a real batch -- it costs one double per output per
// individual per spec. The tests want it, because the per-individual values are what the analytic
// oracles are written against and a plain sum of them must agree with the chunked total.
[[cpp11::register]]
cpp11::list cpp_veil_run(
    cpp11::list specs,
    cpp11::list columns,
    int time_scale,
    bool keep_contributions,
    // Spelled at every call site for the reason given on cpp_veil_aev: cpp11 does not carry a C++
    // default into the R wrapper it generates.
    int threads)
{
    if (threads < 0) { cpp11::stop("`threads` cannot be negative."); }

    const cpp11::strings columnNames = columnNamesOf(columns);
    std::vector<veil::TypeFull> columnTypes;
    std::vector<std::optional<veil::TypeWithConstraints>> constraints;
    veil::ColumnSet set;
    const TextEncoding encoding = prepareColumns(columns, columnNames, set, columnTypes, constraints);

    veil::ExposureColumns exposure;
    exposure.start = exposureColumn(columnNames, columnTypes, "E2R_start", veil::Type::Datey);
    exposure.end = exposureColumn(columnNames, columnTypes, "E2R_end", veil::Type::Datey);
    exposure.died = exposureColumn(columnNames, columnTypes, "E2R_died", veil::Type::Bool);
    exposure.deltaTClicks = veil::validateTimeScaleClicks(time_scale);

    std::vector<AevBuild> builds;
    std::vector<double> dispersions;
    builds.reserve(static_cast<size_t>(specs.size()));
    dispersions.reserve(static_cast<size_t>(specs.size()));
    for (R_xlen_t i = 0; i < specs.size(); ++i)
    {
        const cpp11::list spec = cpp11::as_cpp<cpp11::list>(specs[i]);
        const SEXP mortality = elementOrNull(spec, "mortality");
        if (mortality == R_NilValue)
        {
            cpp11::stop("Specification %d has no `mortality`.", static_cast<int>(i + 1));
        }

        const SEXP overdispersion = elementOrNull(spec, "overdispersion");
        if (overdispersion == R_NilValue)
        {
            cpp11::stop("Specification %d has no `overdispersion`.", static_cast<int>(i + 1));
        }
        dispersions.push_back(checkedOverdispersion(cpp11::as_cpp<double>(overdispersion)));

        builds.push_back(buildAev(
            cpp11::as_cpp<cpp11::list>(mortality),
            elementOrNull(spec, "weight"),
            elementOrNull(spec, "val_similarity"),
            elementOrNull(spec, "val_distance"),
            elementOrNull(spec, "include"),
            columnNames,
            columnTypes,
            constraints,
            encoding.mapping,
            exposure));
    }

    // Taken only once every block exists, because `builds` reallocates as it grows and a pointer into
    // it taken earlier would dangle.
    std::vector<const veil::Block*> blocks;
    blocks.reserve(builds.size());
    for (const AevBuild& build : builds) { blocks.push_back(&build.block); }

    const R_xlen_t records = columns.size() == 0 ? 0 : Rf_xlength(VECTOR_ELT(SEXP(columns), 0));
    const std::vector<const veil::ColumnView*> views = viewsByColumnId(set, columnNames);

    const size_t threadCount = static_cast<size_t>(threads);
    const std::vector<veil::CalculationResult> calculations = veil::runCalculations(
        blocks,
        views,
        static_cast<size_t>(records),
        keep_contributions,
        threadCount,
        userInterruptIsPending);

    // An interrupt stops the whole batch, so one result reporting it means none of them is safe to
    // read. Raised here rather than delivered by R, because the poll consumed the pending interrupt.
    for (const veil::CalculationResult& calculation : calculations)
    {
        if (calculation.interrupted) { cpp11::stop("The veil calculation was interrupted."); }
    }

    const R_xlen_t contributionRows = keep_contributions ? records : 0;

    cpp11::writable::list results;
    for (size_t i = 0; i < calculations.size(); ++i)
    {
        const veil::CalculationResult& calculation = calculations[i];
        const AevBuild& build = builds[i];
        const double dispersion = dispersions[i];

        cpp11::writable::doubles a(contributionRows);
        cpp11::writable::doubles e(contributionRows);
        cpp11::writable::doubles v(contributionRows);
        for (R_xlen_t record = 0; record < contributionRows; ++record)
        {
            const size_t base = static_cast<size_t>(record) * 3;
            a[record] = calculation.contributions[base];
            e[record] = calculation.contributions[base + 1];
            v[record] = varianceFromSecondMoment(calculation.contributions[base + 2], dispersion);
        }

        cpp11::writable::list entry({
            "A"_nm = calculation.totals[0],
            "E"_nm = calculation.totals[1],
            "V"_nm = varianceFromSecondMoment(calculation.totals[2], dispersion),
            "contributions"_nm = cpp11::writable::list({"A"_nm = a, "E"_nm = e, "V"_nm = v}),
            "records_included"_nm = static_cast<int>(calculation.recordsIncluded),
            "slot_evaluations"_nm = static_cast<int>(calculation.slotEvaluations),
            "instruction_count"_nm = static_cast<int>(build.block.body().size()),
            "shared_nodes"_nm = static_cast<int>(build.shared),
            "vector_operand_count"_nm = static_cast<int>(build.layout.vectorOperands),
            "buffer_count"_nm = static_cast<int>(build.layout.buffers),
            "death_only_count"_nm = static_cast<int>(build.layout.deathOnly),
            "dead_instruction_count"_nm = static_cast<int>(build.layout.dead),
        });
        results.push_back(entry);
    }

    return cpp11::writable::list({
        "results"_nm = results,
        "spec_count"_nm = static_cast<int>(builds.size()),

        // Every block divides the same dataset, so one chunk count describes the whole batch. The
        // task count is what threading actually sees, and it is the product.
        "chunk_count"_nm = static_cast<int>(veil::chunkCount(static_cast<size_t>(records))),
        "task_count"_nm =
            static_cast<int>(builds.size() * veil::chunkCount(static_cast<size_t>(records))),

        // THE ONLY OBSERVABLE THAT A THREAD COUNT WAS HONOURED -- see cpp_veil_aev.
        "threads_used"_nm = static_cast<int>(veil::resolveThreadCount(threadCount)),
    });
}

// The same pipeline, additionally reporting the interval worked out for every node so the
// propagation pass can be inspected from R. Intervals are in each node's OWN units, so a datey or
// durationy node reports CLICKS -- see the note at the top of passPropagateIntervals.hpp.
[[cpp11::register]]
cpp11::list cpp_veil_intervals(cpp11::list node, cpp11::list columns)
{
    const cpp11::strings columnNames = columnNamesOf(columns);
    std::vector<veil::TypeFull> columnTypes;
    std::vector<std::optional<veil::TypeWithConstraints>> constraints;
    veil::ColumnSet set;
    const TextEncoding encoding = prepareColumns(columns, columnNames, set, columnTypes, constraints);

    veil::Tree tree;
    veil::ObjStore objs;
    tree.setRoot(ingest(tree, objs, node, columnNames));
    const PipelineResult result =
        runTreePipeline(tree, objs, columnTypes, constraints, encoding.mapping);

    return dumpTree(tree, result.timeVarying, encoding.mapping, &result.intervals);
}
