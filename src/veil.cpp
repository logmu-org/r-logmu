// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

// The veil core translation unit.
//
// veil is plain C++20 with no R dependency: no cpp11, no SEXP, no R runtime calls. That is what
// lets a Python or C# front end reuse this code rather than reimplement it, so nothing R-specific
// may be included from here. The R binding lives outside the core, alongside vec_ops_for_R.cpp.
//
// The frame is header-only, so without this file nothing under src/veil/ would be compiled at all
// and the headers would rot untested. Including them here puts them in the build; the integrity of
// the op table is asserted at compile time inside Op.hpp.

// The data structures.
#include "veil/Block.hpp"
#include "veil/ColumnScan.hpp"
#include "veil/ColumnSet.hpp"
#include "veil/ColumnView.hpp"
#include "veil/DateyInterval.hpp"
#include "veil/Engine.hpp"
#include "veil/Instruction.hpp"
#include "veil/Interpreter.hpp"
#include "veil/MortalityTable.hpp"
#include "veil/Node.hpp"
#include "veil/ObjStore.hpp"
#include "veil/Op.hpp"
#include "veil/Operand.hpp"
#include "veil/RecordChunk.hpp"
#include "veil/ScalarValue.hpp"
#include "veil/StringMapping.hpp"
#include "veil/ThreadPool.hpp"
#include "veil/TimeGrid.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"
#include "veil/TypeSpecificConstraint.hpp"
#include "veil/double_properties.hpp"

// The tree passes. Each header is named for the single pass function it exports, so the two always
// match and they sort together under the `pass` prefix.
#include "veil/passAnnotateTypes.hpp"
#include "veil/passCheckTimeInvariance.hpp"
#include "veil/passFoldConstants.hpp"
#include "veil/passFoldIndicatorSquares.hpp"
#include "veil/passFoldIntervalComparisons.hpp"
#include "veil/passHoistFromIntegrate.hpp"
#include "veil/passInsertCoercions.hpp"
#include "veil/passNarrowComparisons.hpp"
#include "veil/passPropagateIntervals.hpp"
#include "veil/passShareCommonSubtrees.hpp"
#include "veil/passTagTimeVarying.hpp"

// Lowering, which turns the rewritten tree into the executable three-address form.
#include "veil/passLowerToBlock.hpp"

// The TAC passes, which work on the lowered block rather than the tree: what it computes is settled
// by the time these run, and they decide how it is laid out.
#include "veil/passAllocateVectorBuffers.hpp"
#include "veil/passComputeLiveness.hpp"
#include "veil/passComputeSlotDemand.hpp"
#include "veil/passFindDeadInstructions.hpp"

// Recipes: the assemblies that turn a mortality and a weight into the roots of a calculation.
#include "veil/AevRecipe.hpp"
