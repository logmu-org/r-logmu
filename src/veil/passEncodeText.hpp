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
#include "veil/Op.hpp"
#include "veil/StringMapping.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// Turns the strings a front end wrote into the indices the core compares.
//
// A text literal crosses as a string and reaches here still one: `.i$sex == "male"` is a category
// column against text until this pass resolves "male" against the crossing's StringMapping. After it,
// no Type::Text remains anywhere in the tree, which is the postcondition the lowering guard leans on.
//
// WHY THE CORE RESOLVES THE STRING RATHER THAN THE FRONT END. Resolving is one lookup against a table
// of levels; interning is a pass over every row. R does the second, as `factor()`, because it already
// has it in C. Doing the first here as well would have meant a front end rewriting its own syntax
// tree, and would have limited it to a literal written next to its column -- where this reads the
// annotated types and so finds a literal wherever a category value is expected.
//
// AN ABSENT STRING IS ORDINARY, NOT AN ERROR. A typo, or a level in one dataset of a batch and not
// another, gives a comparison that no record can satisfy, so it folds to a constant here and the
// engine never runs it. That is the whole reason `indexOfString` reports absence rather than
// answering with a sentinel index: an index that names no string would then have to be compared
// against, and would collide with the marker for a record that holds no value at all.
//
// RUNS DIRECTLY AFTER passAnnotateTypes, which is what makes the operand types available, and before
// passInsertCoercions, which has nothing to say about a category. Folding a comparison to a bool
// leaves the annotation it already carried correct, so nothing needs re-annotating.

namespace detail
{

inline bool isTextLiteral(const Node& node) noexcept
{
    return isLit(node) && node.type.has_value() && node.type->type == Type::Text
        && std::holds_alternative<std::string>(std::get<LitPayload>(node.payload).value);
}

inline bool isCategory(const Node& node) noexcept
{
    return node.type.has_value() && node.type->type == Type::Category;
}

} // namespace detail

inline void passEncodeText(Tree& tree, const StringMapping& mapping)
{
    // Literals left stranded by a fold: the comparison above them became a constant, so they are
    // still text but no longer reachable. Tracked rather than rewritten, because a literal may be
    // shared and rewriting it would be reaching into a use this pass is not looking at.
    std::vector<char> stranded(tree.size(), 0);

    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        // Re-read rather than held: the loop rewrites nodes, and a reference taken before a rewrite
        // would be describing the payload the node no longer has.
        if (!isCall(tree.at(id))) { continue; }

        const CallPayload call = std::get<CallPayload>(tree.at(id).payload);
        if ((call.op != Op::Eq && call.op != Op::Ne) || call.args.size() != 2) { continue; }

        // Which side is the string and which the category. Ordering comparisons never reach here --
        // passAnnotateTypes refuses them on type -- so equality is the whole of it.
        const NodeId first = call.args[0];
        const NodeId second = call.args[1];

        NodeId literalId = invalidNodeId;
        if (detail::isTextLiteral(tree.at(first)) && detail::isCategory(tree.at(second)))
        {
            literalId = first;
        }
        else if (detail::isTextLiteral(tree.at(second)) && detail::isCategory(tree.at(first)))
        {
            literalId = second;
        }
        else
        {
            continue;
        }

        const std::string& text = std::get<std::string>(std::get<LitPayload>(tree.at(literalId).payload).value);
        const std::optional<int32_t> index = mapping.indexOfString(text);

        if (!index.has_value())
        {
            // No record holds this string, so the answer is settled for every record: never equal,
            // and therefore always unequal.
            const bool answer = call.op == Op::Ne;
            tree.at(id).payload.emplace<LitPayload>(LitPayload{LiteralValue(answer), TypeFull::createBool()});
            stranded[literalId] = 1;
            continue;
        }

        // Encoded on the literal itself, so a literal reached by two comparisons is resolved once and
        // read twice: the second visit finds a category rather than text and passes over it. `%in%`
        // desugars to one comparison per element, so more than one use is the ordinary shape here.
        //
        // `emplace` rather than assignment on both: TypeFull holds its members const, which makes it
        // constructible but not assignable.
        const TypeFull codeType = TypeFull::createCategory(static_cast<int>(mapping.size()));
        tree.at(literalId).payload.emplace<LitPayload>(LitPayload{LiteralValue(static_cast<int>(*index)), codeType});
        tree.at(literalId).type.emplace(codeType);
    }

    // NOTHING MAY STILL BE TEXT. A text literal that met no category value, or a character column
    // that never became a factor, would otherwise reach lowering and be refused there on storage --
    // a complaint about how veil holds a value, to a user who asked a question about strings.
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        const Node& node = tree.at(id);
        if (stranded[id] != 0) { continue; }
        if (!node.type.has_value() || node.type->type != Type::Text) { continue; }

        if (isField(node))
        {
            throw std::runtime_error(
                "veil: column `" + std::get<FieldPayload>(node.payload).name + "` holds text that has "
                "not been given category levels; a text column must arrive as a factor.");
        }
        throw std::runtime_error(
            "veil: a text value can only be compared with a category value, using `==` or `!=`.");
    }
}

} // namespace veil
