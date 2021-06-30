// Author: Enrico Guiraud, Danilo Piparo CERN  12/2016

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TDFUTILS
#define ROOT_TDFUTILS

#include "ROOT/TDataSource.hxx" // ColumnName2ColumnTypeName
#include "ROOT/TypeTraits.hxx"
#include "ROOT/TArrayBranch.hxx"
#include "Compression.h"
#include "TH1.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <array>
#include <cstddef> // std::size_t
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <type_traits> // std::decay
#include <vector>

class TTree;
class TTreeReader;

namespace ROOT {
namespace Experimental {
namespace TDF {
/// A collection of options to steer the creation of the dataset on file
struct TSnapshotOptions {
   using ECAlgo = ::ROOT::ECompressionAlgorithm;
   TSnapshotOptions() = default;
   TSnapshotOptions(const TSnapshotOptions &) = default;
   TSnapshotOptions(TSnapshotOptions &&) = default;
   TSnapshotOptions(std::string_view mode, ECAlgo comprAlgo, int comprLevel, int autoFlush, int splitLevel)
      : fMode(mode), fCompressionAlgorithm(comprAlgo), fCompressionLevel{comprLevel}, fAutoFlush(autoFlush),
        fSplitLevel(splitLevel)
   {
   }
   std::string fMode = "RECREATE";             //< Mode of creation of output file
   ECAlgo fCompressionAlgorithm = ROOT::kZLIB; //< Compression algorithm of output file
   int fCompressionLevel = 1;                  //< Compression level of output file
   int fAutoFlush = 0;                         //< AutoFlush value for output tree
   int fSplitLevel = 99;                       //< Split level of output tree
};
}
}
}

/// \cond HIDDEN_SYMBOLS

namespace ROOT {

// fwd declaration for IsV7Hist
namespace Experimental {
template <int D, typename P, template <int, typename, template <typename> class> class... S>
class THist;

} // ns Experimental

namespace Detail {
namespace TDF {
using ColumnNames_t = std::vector<std::string>;

// fwd decl for ColumnName2ColumnTypeName
class TCustomColumnBase;
// fwd decl for FindUnknownColumns
class TLoopManager;

// type used for tag dispatching
struct TInferType {
};

} // end ns Detail
} // end ns TDF

namespace Internal {
namespace TDF {
using namespace ROOT::TypeTraits;
using namespace ROOT::Detail::TDF;
using namespace ROOT::Experimental::TDF;

class TIgnoreErrorLevelRAII {
private:
   int fCurIgnoreErrorLevel = gErrorIgnoreLevel;

public:
   TIgnoreErrorLevelRAII(int errorIgnoreLevel);
   ~TIgnoreErrorLevelRAII();
};

/// Compile-time integer sequence generator
/// e.g. calling GenStaticSeq<3>::type() instantiates a StaticSeq<0,1,2>
template <int...>
struct StaticSeq {
};

template <int N, int... S>
struct GenStaticSeq : GenStaticSeq<N - 1, N - 1, S...> {
};

template <int... S>
struct GenStaticSeq<0, S...> {
   using type = StaticSeq<S...>;
};

template <int... S>
using GenStaticSeq_t = typename GenStaticSeq<S...>::type;

// return wrapper around f that prepends an `unsigned int slot` parameter
template <typename R, typename F, typename... Args>
std::function<R(unsigned int, Args...)> AddSlotParameter(F &f, TypeList<Args...>)
{
   return [f](unsigned int, Args... a) -> R { return f(a...); };
}

template <typename BranchType, typename... Rest>
struct TNeedJitting {
   static constexpr bool value = TNeedJitting<Rest...>::value;
};

template <typename... Rest>
struct TNeedJitting<TInferType, Rest...> {
   static constexpr bool value = true;
};

template <typename T>
struct TNeedJitting<T> {
   static constexpr bool value = false;
};

template <>
struct TNeedJitting<TInferType> {
   static constexpr bool value = true;
};

using TVBPtr_t = std::shared_ptr<TTreeReaderValueBase>;
using TVBVec_t = std::vector<TVBPtr_t>;

const std::type_info &TypeName2TypeID(const std::string &name);

std::string TypeID2TypeName(const std::type_info &id);

std::string
ColumnName2ColumnTypeName(const std::string &colName, TTree *, TCustomColumnBase *, TDataSource * = nullptr);

char TypeName2ROOTTypeName(const std::string &b);

const char *ToConstCharPtr(const char *s);
const char *ToConstCharPtr(const std::string &s);
unsigned int GetNSlots();

/// Choose between TTreeReader{Array,Value} depending on whether the branch type
/// T is a `TArrayBranch<T>` or any other type (respectively).
template <typename T>
struct TReaderValueOrArray {
   using Proxy_t = TTreeReaderValue<T>;
};

template <typename T>
struct TReaderValueOrArray<TArrayBranch<T>> {
   using Proxy_t = TTreeReaderArray<T>;
};

template <typename T>
using ReaderValueOrArray_t = typename TReaderValueOrArray<T>::Proxy_t;

/// Initialize a tuple of TColumnValues.
/// For real TTree branches a TTreeReader{Array,Value} is built and passed to the
/// TColumnValue. For temporary columns a pointer to the corresponding variable
/// is passed instead.
template <typename TDFValueTuple, int... S>
void InitTDFValues(unsigned int slot, TDFValueTuple &valueTuple, TTreeReader *r, const ColumnNames_t &bn,
                   const ColumnNames_t &tmpbn,
                   const std::map<std::string, std::shared_ptr<TCustomColumnBase>> &customCols, StaticSeq<S...>)
{
   // isTmpBranch has length bn.size(). Elements are true if the corresponding
   // branch is a temporary branch created with Define, false if they are
   // actual branches present in the TTree.
   std::array<bool, sizeof...(S)> isTmpColumn;
   for (auto i = 0u; i < isTmpColumn.size(); ++i)
      isTmpColumn[i] = std::find(tmpbn.begin(), tmpbn.end(), bn.at(i)) != tmpbn.end();

   // hack to expand a parameter pack without c++17 fold expressions.
   // The statement defines a variable with type std::initializer_list<int>, containing all zeroes, and SetTmpColumn or
   // SetProxy are conditionally executed as the braced init list is expanded. The final ... expands S.
   std::initializer_list<int> expander{(isTmpColumn[S]
                                           ? std::get<S>(valueTuple).SetTmpColumn(slot, customCols.at(bn.at(S)).get())
                                           : std::get<S>(valueTuple).MakeProxy(r, bn.at(S)),
                                        0)...};
   (void)expander; // avoid "unused variable" warnings for expander on gcc4.9
   (void)slot;     // avoid _bogus_ "unused variable" warnings for slot on gcc 4.9
   (void)r;        // avoid "unused variable" warnings for r on gcc5.2
}

template <typename Filter>
void CheckFilter(Filter &)
{
   using FilterRet_t = typename TDF::CallableTraits<Filter>::ret_type;
   static_assert(std::is_same<FilterRet_t, bool>::value, "filter functions must return a bool");
}

ColumnNames_t GetBranchNames(TTree &t);

void CheckCustomColumn(std::string_view definedCol, TTree *treePtr, const ColumnNames_t &customCols,
                       const ColumnNames_t &dataSourceColumns);

///////////////////////////////////////////////////////////////////////////////
/// Check that the callable passed to TInterface::Reduce:
/// - takes exactly two arguments of the same type
/// - has a return value of the same type as the arguments
template <typename F, typename T>
void CheckReduce(F &, TypeList<T, T>)
{
   using ret_type = typename CallableTraits<F>::ret_type;
   static_assert(std::is_same<ret_type, T>::value, "reduce function must have return type equal to argument type");
   return;
}

///////////////////////////////////////////////////////////////////////////////
/// This overload of CheckReduce is called if T is not a TypeList<T,T>
template <typename F, typename T>
void CheckReduce(F &, T)
{
   static_assert(sizeof(F) == 0, "reduce function must take exactly two arguments of the same type");
}

///////////////////////////////////////////////////////////////////////////////
/// Check as many template parameters were passed as the number of column names, throw if this is not the case.
void CheckSnapshot(unsigned int nTemplateParams, unsigned int nColumnNames);

/// Return local BranchNames or default BranchNames according to which one should be used
const ColumnNames_t SelectColumns(unsigned int nArgs, const ColumnNames_t &bl, const ColumnNames_t &defBl);

/// Check whether column names refer to a valid branch of a TTree or have been `Define`d. Return invalid column names.
ColumnNames_t FindUnknownColumns(const ColumnNames_t &requiredCols, TTree *tree, const ColumnNames_t &definedCols,
                                 const ColumnNames_t &dataSourceColumns);

// clang-format off
namespace ActionTypes {
struct Histo1D {};
struct Histo2D {};
struct Histo3D {};
struct Profile1D {};
struct Profile2D {};
struct Min {};
struct Max {};
struct Sum {};
struct Mean {};
struct Fill {};
}
// clang-format on

/// Check whether a histogram type is a classic or v7 histogram.
template <typename T>
struct IsV7Hist : public std::false_type {
   static_assert(std::is_base_of<TH1, T>::value, "not implemented for this type");
};

template <int D, typename P, template <int, typename, template <typename> class> class... S>
struct IsV7Hist<ROOT::Experimental::THist<D, P, S...>> : public std::true_type {
};

template <typename T, bool ISV7HISTO = IsV7Hist<T>::value>
struct HistoUtils {
   static void SetCanExtendAllAxes(T &h) { h.SetCanExtend(::TH1::kAllAxes); }
   static bool HasAxisLimits(T &h)
   {
      auto xaxis = h.GetXaxis();
      return !(xaxis->GetXmin() == 0. && xaxis->GetXmax() == 0.);
   }
};

template <typename T>
struct HistoUtils<T, true> {
   static void SetCanExtendAllAxes(T &) {}
   static bool HasAxisLimits(T &) { return true; }
};

/// `type` is TypeList if MustRemove is false, otherwise it is a TypeList with the first type removed
template <bool MustRemove, typename TypeList>
struct RemoveFirstParameterIf {
   using type = TypeList;
};

template <typename TypeList>
struct RemoveFirstParameterIf<true, TypeList> {
   using type = RemoveFirstParameter_t<TypeList>;
};

template <bool MustRemove, typename TypeList>
struct RemoveFirstTwoParametersIf {
   using type = TypeList;
};

template <typename TypeList>
struct RemoveFirstTwoParametersIf<true, TypeList> {
   using typeTmp = typename RemoveFirstParameterIf<true, TypeList>::type;
   using type = typename RemoveFirstParameterIf<true, typeTmp>::type;
};

bool IsInternalColumn(std::string_view colName);

// Check if a condition is true for all types
template <bool...>
struct TBoolPack;
template <bool... bs>
using IsTrueForAllImpl_t = typename std::is_same<TBoolPack<bs..., true>, TBoolPack<true, bs...>>;

template <bool... Conditions>
struct TEvalAnd {
   static constexpr bool value = IsTrueForAllImpl_t<Conditions...>::value;
};

// Check if a class is a specialisation of stl containers templates

template <typename>
struct IsVector_t : std::false_type {
};

template <typename T>
struct IsVector_t<std::vector<T>> : std::true_type {
};

template <typename>
struct IsList_t : std::false_type {
};

template <typename T>
struct IsList_t<std::list<T>> : std::true_type {
};

template <typename>
struct IsDeque_t : std::false_type {
};

template <typename T>
struct IsDeque_t<std::deque<T>> : std::true_type {
};

template <typename>
struct IsTArrayBranch_t : std::false_type {
};

template <typename T>
struct IsTArrayBranch_t<ROOT::Experimental::TDF::TArrayBranch<T>> : std::true_type {
};

// Check the value_type type of a type with a SFINAE to allow compilation in presence
// fundamental types
template <typename T, bool IsContainer = IsContainer<typename std::decay<T>::type>::value>
struct ValueType {
   using value_type = typename T::value_type;
};

template <typename T>
struct ValueType<T, false> {
   using value_type = T;
};

template <typename T>
struct ValueType<ROOT::Experimental::TDF::TArrayBranch<T>, false> {
   using value_type = T;
};

} // end NS TDF
} // end NS Internal
} // end NS ROOT

/// \endcond

#endif // TDFUTILS
