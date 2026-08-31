// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>
#include "veil/Node.hpp"
#include "veil/ObjStore.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// The type annotation pass. It fills every node's base type and nothing else -- no structural
// change. It follows the datey package's operator rules (https://r-datey.logmu.org/reference/ops.html):
// a number anywhere in a value operator makes the result double (datey / durationy read as years),
// and pure datey / durationy arithmetic uses the exact integer algebra.
//
// Kept deliberately non-mutating so it stays a re-runnable "type everything" utility that a later
// pass can call again after it creates nodes. Coercion-node insertion is a separate pass that reads
// the types written here.

// A readable name for a type, for diagnostic messages.
inline std::string typeLabel(const TypeFull& type)
{
    switch (type.type)
    {
        case Type::Bool: return "bool";
        case Type::Double: return "double";
        case Type::Datey: return "datey";
        case Type::DateyInterval: return "datey_interval";
        case Type::Durationy: return "durationy";
        case Type::Text: return "text";
        case Type::Category: return "category";
    }
    return "unknown";
}

namespace detail
{

// The double-valued family: a plain number, or a bool that reads as 0 / 1.
inline bool isDoubleFamily(Type t) noexcept { return t == Type::Double || t == Type::Bool; }

// The integer clicks on the annual grid.
inline bool isClickInt(Type t) noexcept { return t == Type::Datey || t == Type::Durationy; }

// Everything that can take part in arithmetic as a number.
inline bool isNumeric(Type t) noexcept { return isDoubleFamily(t) || isClickInt(t); }

// The two forms one category value takes: the string a front end wrote, and the index the core
// compares. `passEncodeText` turns the first into the second, and runs straight after this pass.
inline bool isTextual(Type t) noexcept { return t == Type::Text || t == Type::Category; }

// Two operands are comparable when they are the same type, or one is a plain number that widens the
// other (datey / durationy to years, bool to 0 / 1).
//
// TEXT AGAINST CATEGORY IS THE UNENCODED FORM OF CATEGORY AGAINST CATEGORY. `.i$sex == "male"` is a
// column of indices against a string, and stays that way until `passEncodeText` resolves the string
// against the crossing's StringMapping. Admitting it here rather than after encoding is what lets the
// encoding pass read the annotated types it needs to know which literals are category values at all.
inline bool comparable(Type a, Type b) noexcept
{
    if (a == b) { return true; }
    if (a == Type::Double && isNumeric(b)) { return true; }
    if (b == Type::Double && isNumeric(a)) { return true; }
    if (isTextual(a) && isTextual(b)) { return true; }
    return false;
}

inline std::string monikerString(Op op) { return std::string(moniker(op)); }

[[noreturn]] inline void undefinedFor(Op op, const std::vector<TypeFull>& args)
{
    std::string message = "veil: `" + monikerString(op) + "` is not defined for ";
    for (size_t i = 0; i < args.size(); ++i)
    {
        if (i != 0) { message += (i + 1 == args.size()) ? " and " : ", "; }
        message += typeLabel(args[i]);
    }
    message += ".";
    throw std::runtime_error(message);
}

// The exact datey / durationy algebra, reached only when every operand is a click integer.
inline TypeFull clickArithmetic(Op op, const std::vector<TypeFull>& args)
{
    const Type a = args[0].type;
    switch (op)
    {
        case Op::Pos:
            return args[0];
        case Op::Neg:
        case Op::Abs:
            if (a == Type::Durationy) { return TypeFull::createDurationy(); }
            break; // Negating or taking the magnitude of a date is undefined.
        case Op::Add:
        {
            const Type b = args[1].type;
            if (a == Type::Durationy && b == Type::Durationy) { return TypeFull::createDurationy(); }
            if ((a == Type::Datey && b == Type::Durationy) || (a == Type::Durationy && b == Type::Datey))
            {
                return TypeFull::createDatey();
            }
            break; // datey + datey is an error.
        }
        case Op::Sub:
        {
            const Type b = args[1].type;
            if (a == Type::Datey && b == Type::Datey) { return TypeFull::createDurationy(); }
            if (a == Type::Datey && b == Type::Durationy) { return TypeFull::createDatey(); }
            if (a == Type::Durationy && b == Type::Durationy) { return TypeFull::createDurationy(); }
            break; // durationy - datey is an error.
        }
        case Op::Min:
        case Op::Max:
            if (a == args[1].type) { return args[0]; }
            break; // A date and a duration have no common extreme.
        case Op::Clamp:
            if (a == args[1].type && a == args[2].type) { return args[0]; }
            break;
        default:
            break; // mul, idiv, imod and the rest need a number.
    }
    undefinedFor(op, args);
}

// The multiplicative and additive operators, following datey: a number anywhere yields double, else
// the pure click algebra.
inline TypeFull arithmeticResult(Op op, const std::vector<TypeFull>& args)
{
    for (const TypeFull& arg : args)
    {
        if (!isNumeric(arg.type))
        {
            throw std::runtime_error("veil: `" + monikerString(op) + "` cannot be applied to "
                                     + typeLabel(arg) + ".");
        }
    }

    for (const TypeFull& arg : args)
    {
        if (isDoubleFamily(arg.type)) { return TypeFull::createDouble(); }
    }

    return clickArithmetic(op, args);
}

// The transcendental and always-double operators: every operand must be numeric; the result is
// double. datey / durationy operands read as years.
inline TypeFull doubleResult(Op op, const std::vector<TypeFull>& args)
{
    for (const TypeFull& arg : args)
    {
        if (!isNumeric(arg.type))
        {
            throw std::runtime_error("veil: `" + monikerString(op) + "` expects numeric operands, got "
                                     + typeLabel(arg) + ".");
        }
    }
    return TypeFull::createDouble();
}

inline TypeFull boolResult(Op op, const std::vector<TypeFull>& args)
{
    if (opInfo(op).category == OpCategory::Logical)
    {
        for (const TypeFull& arg : args)
        {
            if (arg.type != Type::Bool)
            {
                throw std::runtime_error("veil: `" + monikerString(op) + "` expects logical operands, got "
                                         + typeLabel(arg) + ".");
            }
        }
        return TypeFull::createBool();
    }

    if (op == Op::IsNa)
    {
        const Type t = args[0].type;
        if (t != Type::Double && !isClickInt(t))
        {
            throw std::runtime_error("veil: `is_na` applies to double, datey or durationy, not "
                                     + typeLabel(args[0]) + ".");
        }
        return TypeFull::createBool();
    }

    // TEXT HAS NO ORDERING, so `<` on it is refused here on type -- `comparable` says only that two
    // operands may be compared, not by which operator. Refused at annotation rather than at lowering
    // because it is a fact about the type.
    //
    // CATEGORY IS BANNED ALONGSIDE TEXT AND IS THE CASE THAT ACTUALLY FIRES. A text column arrives as
    // a factor, so `.i$sex < "male"` is category against text and never text against text; leaving
    // this on Text alone would have let every real expression through to lowering, which refuses it
    // with a message about storage rather than about meaning.
    //
    // THE MESSAGE SAYS TEXT WHATEVER THE OPERAND TYPES SAY (Tim, 2026-08-11). A user writing this
    // wrote strings and has never heard the word category: it is the name of the encoded form they
    // never see. Naming the type they wrote is what makes the complaint actionable.
    //
    // CHECKED BEFORE `comparable`, so that the answer to `.i$sex < .i$age` is that text cannot be
    // ordered rather than that the two types do not compare -- the ordering is the user's mistake and
    // the mismatch merely follows from it.
    if (isOrderingComparison(op) && (isTextual(args[0].type) || isTextual(args[1].type)))
    {
        throw std::runtime_error("veil: `" + monikerString(op) + "` cannot order text; text records "
                                 "which value, not where it sits in a sequence, so only `==` and "
                                 "`!=` apply.");
    }

    // A comparison.
    if (!comparable(args[0].type, args[1].type))
    {
        throw std::runtime_error("veil: `" + monikerString(op) + "` cannot compare " + typeLabel(args[0])
                                 + " with " + typeLabel(args[1]) + ".");
    }
    return TypeFull::createBool();
}

// The common type of a select's two branches. Equal types agree; otherwise a branch widens only when
// a plain number meets a number.
// A SELECT MAY PRODUCE A CATEGORY BUT MAY NOT CHOOSE BETWEEN WRITTEN STRINGS (ruled 2026-08-11).
//
// `ifelse(.i$smoker, .i$sex, .i$spouse)` is allowed and works: both branches are columns, so both are
// already indices into the crossing's StringMapping and the select is an ordinary choice between two
// values of one type. That shape is a real one -- test against the member's own sex or their
// spouse's, depending on a flag.
//
// `ifelse(.i$smoker, "heavy", "light")` is refused. The mapping is built from the COLUMNS, at the
// boundary, before any expression is typed -- `Category.max` on every column's type is its size -- so
// a string that no column holds has no index and cannot be given one without rebuilding the mapping
// after ingest and re-typing everything against it. That is real machinery.
//
// What it would buy is nothing. A category value's only operation is equality, so every use of a
// computed one is a comparison, and every such comparison is the plain disjunction written out:
// `ifelse(c, "a", "b") == .i$x` is `(c & .i$x == "a") | (!c & .i$x == "b")`. No expressive power is
// lost, and the rule stays statable in one line.
//
// TO REVERSE THIS: give the mapping the tree's text literals as well as the columns' levels, which
// means finding them before `rColumnType` runs. Nothing else here would need to change.
inline TypeFull branchCommonType(const TypeFull& x, const TypeFull& y)
{
    if (x.type == Type::Text || y.type == Type::Text)
    {
        throw std::runtime_error("veil: `select` cannot choose between text values; compare each "
                                 "branch instead, as `(c & x == \"a\") | (!c & x == \"b\")`.");
    }

    if (x == y) { return x; }
    if (isDoubleFamily(x.type) && isDoubleFamily(y.type)) { return TypeFull::createDouble(); }
    if (x.type == Type::Double && isClickInt(y.type)) { return TypeFull::createDouble(); }
    if (y.type == Type::Double && isClickInt(x.type)) { return TypeFull::createDouble(); }
    throw std::runtime_error("veil: `select` branches have incompatible types " + typeLabel(x) + " and "
                             + typeLabel(y) + ".");
}

inline TypeFull selectResult(const std::vector<TypeFull>& args)
{
    if (args[0].type != Type::Bool)
    {
        throw std::runtime_error("veil: `select` needs a logical condition, got " + typeLabel(args[0]) + ".");
    }
    return branchCommonType(args[1], args[2]);
}

// An obj node's type. An include resolves to a datey_interval; a mortality table is left untyped,
// because it is only ever the first argument of vector_log_mu, which is checked directly.
inline std::optional<TypeFull> objType(const ObjPayload& payload, const ObjStore& objs)
{
    const Obj& obj = objs.at(payload.obj);
    if (std::holds_alternative<Include>(obj)) { return TypeFull::createDateyInterval(); }
    if (std::holds_alternative<MortalityConst>(obj)) { return TypeFull::createDouble(); }
    return std::nullopt; // MortalityTable.
}

// The value type of an operand. Everything except a mortality table has one by the time a parent is
// reached; a table appearing anywhere but vector_log_mu's first argument is caught here.
inline TypeFull operandType(const Tree& tree, NodeId argId)
{
    const std::optional<TypeFull>& type = tree.at(argId).type;
    if (!type.has_value())
    {
        throw std::runtime_error("veil: a mortality table can only be the first argument of `vector_log_mu`.");
    }
    return *type;
}

inline TypeFull callType(const Tree& tree, const CallPayload& call, const ObjStore& objs)
{
    switch (call.op)
    {
        case Op::VectorLogMu:
        {
            const Node& tableNode = tree.at(call.args[0]);
            const bool okTable = isObj(tableNode)
                && std::holds_alternative<MortalityTable>(objs.at(std::get<ObjPayload>(tableNode.payload).obj));
            if (!okTable)
            {
                throw std::runtime_error("veil: `vector_log_mu` expects a mortality table as its first argument.");
            }
            const TypeFull birth = operandType(tree, call.args[1]);
            if (birth.type != Type::Datey)
            {
                throw std::runtime_error("veil: `vector_log_mu` expects a datey `birth`, got " + typeLabel(birth) + ".");
            }
            return TypeFull::createDouble();
        }
        case Op::VectorT:
        case Op::VectorDurn:
        case Op::Broadcast:
            return TypeFull::createDouble(); // A double time vector; operands are plain numerics.
        default:
            break;
    }

    std::vector<TypeFull> args;
    args.reserve(call.args.size());
    for (const NodeId argId : call.args) { args.push_back(operandType(tree, argId)); }

    switch (opInfo(call.op).resultRule)
    {
        case ResultRule::AlwaysBool: return boolResult(call.op, args);
        case ResultRule::AlwaysDouble: return doubleResult(call.op, args);
        case ResultRule::SameAsBranches: return selectResult(args);
        case ResultRule::SameAsArgs: return arithmeticResult(call.op, args);
        case ResultRule::AlwaysVector: return TypeFull::createDouble(); // Handled above; defensive.
    }
    throw std::runtime_error("veil: `" + monikerString(call.op) + "` has no result rule.");
}

} // namespace detail

// Fills the type of every node whose type is not already set, in a single forward pass. The builders
// add a node's children before the node itself, so every argument has a smaller id and is already
// typed by the time its parent is reached.
inline void passAnnotateTypes(Tree& tree, const std::vector<TypeFull>& columnTypes, const ObjStore& objs)
{
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        Node& node = tree.at(id);
        if (node.type.has_value()) { continue; } // Lits are typed by the builders.

        // TypeFull has const members, so optional<TypeFull> is not assignable; emplace constructs
        // the value in place.
        if (isTime(node))
        {
            node.type.emplace(TypeFull::createDatey());
        }
        else if (isField(node))
        {
            const FieldPayload& field = std::get<FieldPayload>(node.payload);
            if (field.column >= columnTypes.size())
            {
                throw std::runtime_error("veil: field `.i$" + field.name + "` refers to an unknown column.");
            }
            node.type.emplace(columnTypes[field.column]);
        }
        else if (isObj(node))
        {
            // A mortality table stays untyped (objType yields no value), which is left as it is.
            if (const std::optional<TypeFull> type = detail::objType(std::get<ObjPayload>(node.payload), objs))
            {
                node.type.emplace(*type);
            }
        }
        else // A call.
        {
            const CallPayload& call = std::get<CallPayload>(node.payload);
            for (const NodeId argId : call.args)
            {
                if (argId >= id) { throw std::logic_error("veil: a call's argument was built after it."); }
            }
            node.type.emplace(detail::callType(tree, call, objs));
        }
    }
}

} // namespace veil
