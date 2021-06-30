// Author: Enrico Guiraud, Danilo Piparo CERN  03/2017

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TDFNODES
#define ROOT_TDFNODES

#include "ROOT/TypeTraits.hxx"
#include "ROOT/TDataSource.hxx"
#include "ROOT/TDFUtils.hxx"
#include "ROOT/TArrayBranch.hxx"
#include "ROOT/TSpinMutex.hxx"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TError.h"

#include <map>
#include <numeric> // std::accumulate (PrintReport), std::iota (TSlotStack)
#include <string>
#include <tuple>
#include <cassert>
#include <climits>
#include <deque> // std::vector substitute in case of vector<bool>
#include <functional>

namespace ROOT {

namespace Internal {
namespace TDF {
class TActionBase;

// This is an helper class to allow to pick a slot without resorting to a map
// indexed by thread ids.
// WARNING: this class does not work as a regular stack. The size is
// fixed at construction time and no blocking is foreseen.
class TSlotStack {
private:
   unsigned int &GetCount()
   {
      TTHREAD_TLS(unsigned int) count = 0U;
      return count;
   }
   unsigned int &GetIndex()
   {
      TTHREAD_TLS(unsigned int) index = UINT_MAX;
      return index;
   }
   unsigned int fCursor;
   std::vector<unsigned int> fBuf;
   ROOT::TSpinMutex fMutex;

public:
   TSlotStack() = delete;
   TSlotStack(unsigned int size) : fCursor(size), fBuf(size) { std::iota(fBuf.begin(), fBuf.end(), 0U); }
   void ReturnSlot(unsigned int slotNumber);
   unsigned int GetSlot();
};
}
}

namespace Detail {
namespace TDF {
using namespace ROOT::TypeTraits;
namespace TDFInternal = ROOT::Internal::TDF;

// forward declarations for TLoopManager
using ActionBasePtr_t = std::shared_ptr<TDFInternal::TActionBase>;
using ActionBaseVec_t = std::vector<ActionBasePtr_t>;
class TCustomColumnBase;
using TCustomColumnBasePtr_t = std::shared_ptr<TCustomColumnBase>;
class TFilterBase;
using FilterBasePtr_t = std::shared_ptr<TFilterBase>;
using FilterBaseVec_t = std::vector<FilterBasePtr_t>;
class TRangeBase;
using RangeBasePtr_t = std::shared_ptr<TRangeBase>;
using RangeBaseVec_t = std::vector<RangeBasePtr_t>;

class TLoopManager : public std::enable_shared_from_this<TLoopManager> {
   using TDataSource = ROOT::Experimental::TDF::TDataSource;
   enum class ELoopType { kROOTFiles, kROOTFilesMT, kNoFiles, kNoFilesMT, kDataSource, kDataSourceMT };

   using Callback_t = std::function<void(unsigned int)>;
   class TCallback {
      const Callback_t fFun;
      const ULong64_t fEveryN;
      std::vector<ULong64_t> fCounters;

   public:
      TCallback(ULong64_t everyN, Callback_t &&f, unsigned int nSlots)
         : fFun(std::move(f)), fEveryN(everyN), fCounters(nSlots, 0ull) {}

      void operator()(unsigned int slot)
      {
         auto &c = fCounters[slot];
         ++c;
         if (c == fEveryN) {
            c = 0ull;
            fFun(slot);
         }
      }
   };

   class TOneTimeCallback {
      const Callback_t fFun;
      std::vector<int> fHasBeenCalled; // std::vector<bool> is thread-unsafe for our purposes (and generally evil)

   public:
      TOneTimeCallback(Callback_t &&f, unsigned int nSlots) : fFun(std::move(f)), fHasBeenCalled(nSlots, 0) {}

      void operator()(unsigned int slot)
      {
         if (fHasBeenCalled[slot] == 1)
            return;
         fFun(slot);
         fHasBeenCalled[slot] = 1;
      }
   };

   ActionBaseVec_t fBookedActions;
   FilterBaseVec_t fBookedFilters;
   FilterBaseVec_t fBookedNamedFilters; ///< Contains a subset of fBookedFilters, i.e. only the named filters
   std::map<std::string, TCustomColumnBasePtr_t> fBookedCustomColumns;
   ColumnNames_t fCustomColumnNames; ///< Contains names of all custom columns defined in the functional graph.
   RangeBaseVec_t fBookedRanges;
   std::vector<std::shared_ptr<bool>> fResProxyReadiness;
   ::TDirectory *const fDirPtr{nullptr};
   std::shared_ptr<TTree> fTree{nullptr}; //< Shared pointer to the input TTree/TChain. It does not own the pointee if
   // the TTree/TChain was passed directly as an argument to TDataFrame's ctor (in
   // which case we let users retain ownership).
   const ColumnNames_t fDefaultColumns;
   const ULong64_t fNEmptyEntries{0};
   const unsigned int fNSlots{1};
   bool fMustRunNamedFilters{true};
   unsigned int fNChildren{0};      ///< Number of nodes of the functional graph hanging from this object
   unsigned int fNStopsReceived{0}; ///< Number of times that a children node signaled to stop processing entries.
   const ELoopType fLoopType; ///< The kind of event loop that is going to be run (e.g. on ROOT files, on no files)
   std::string fToJit;        ///< string containing all `BuildAndBook` actions that should be jitted before running
   const std::unique_ptr<TDataSource> fDataSource; ///< Owning pointer to a data-source object. Null if no data-source
   ColumnNames_t fDefinedDataSourceColumns;        ///< List of data-source columns that have been `Define`d so far
   std::map<std::string, std::string> fAliasColumnNameMap; ///< ColumnNameAlias-columnName pairs
   std::vector<TCallback> fCallbacks; ///< Registered callbacks
   std::vector<TOneTimeCallback> fCallbacksOnce; ///< Registered callbacks to invoke just once before running the loop

   void RunEmptySourceMT();
   void RunEmptySource();
   void RunTreeProcessorMT();
   void RunTreeReader();
   void RunDataSourceMT();
   void RunDataSource();
   void RunAndCheckFilters(unsigned int slot, Long64_t entry);
   void InitNodeSlots(TTreeReader *r, unsigned int slot);
   void InitNodes();
   void CleanUpNodes();
   void CleanUpTask(unsigned int slot);
   void JitActions();
   void EvalChildrenCounts();

public:
   TLoopManager(TTree *tree, const ColumnNames_t &defaultBranches);
   TLoopManager(ULong64_t nEmptyEntries);
   TLoopManager(std::unique_ptr<TDataSource> ds, const ColumnNames_t &defaultBranches);
   TLoopManager(const TLoopManager &) = delete;
   TLoopManager &operator=(const TLoopManager &) = delete;

   void Run();
   TLoopManager *GetImplPtr();
   std::shared_ptr<TLoopManager> GetSharedPtr() { return shared_from_this(); }
   const ColumnNames_t &GetDefaultColumnNames() const;
   const ColumnNames_t &GetCustomColumnNames() const { return fCustomColumnNames; };
   TTree *GetTree() const;
   TCustomColumnBase *GetBookedBranch(const std::string &name) const;
   const std::map<std::string, TCustomColumnBasePtr_t> &GetBookedColumns() const { return fBookedCustomColumns; }
   ::TDirectory *GetDirectory() const;
   ULong64_t GetNEmptyEntries() const { return fNEmptyEntries; }
   TDataSource *GetDataSource() const { return fDataSource.get(); }
   void Book(const ActionBasePtr_t &actionPtr);
   void Book(const FilterBasePtr_t &filterPtr);
   void Book(const TCustomColumnBasePtr_t &branchPtr);
   void Book(const std::shared_ptr<bool> &branchPtr);
   void Book(const RangeBasePtr_t &rangePtr);
   bool CheckFilters(int, unsigned int);
   unsigned int GetNSlots() const { return fNSlots; }
   bool MustRunNamedFilters() const { return fMustRunNamedFilters; }
   void Report() const;
   /// End of recursive chain of calls, does nothing
   void PartialReport() const {}
   void SetTree(const std::shared_ptr<TTree> &tree) { fTree = tree; }
   void IncrChildrenCount() { ++fNChildren; }
   void StopProcessing() { ++fNStopsReceived; }
   void Jit(const std::string &s) { fToJit.append(s); }
   const ColumnNames_t &GetDefinedDataSourceColumns() const { return fDefinedDataSourceColumns; }
   void AddDataSourceColumn(std::string_view name) { fDefinedDataSourceColumns.emplace_back(name); }
   void AddColumnAlias(const std::string &alias, const std::string &colName) { fAliasColumnNameMap[alias] = colName; }
   const std::map<std::string, std::string> &GetAliasMap() const { return fAliasColumnNameMap; }
   void RegisterCallback(ULong64_t everyNEvents, std::function<void(unsigned int)> &&f);
};
} // end ns TDF
} // end ns Detail

namespace Internal {
namespace TDF {
using namespace ROOT::Detail::TDF;

/**
\class ROOT::Internal::TDF::TColumnValue
\ingroup dataframe
\brief Helper class that updates and returns TTree branches as well as TDataFrame temporary columns
\tparam T The type of the column

TDataFrame nodes must access two different types of values during the event loop:
values of real branches, for which TTreeReader{Values,Arrays} act as proxies, or
temporary columns whose values are generated on the fly. While the type of the
value is known at compile time (or just-in-time), it is only at runtime that nodes
can check whether a certain value is generated on the fly or not.

TColumnValue abstracts this difference by providing the same interface for
both cases and handling the reading or generation of new values transparently.
Only one of the two data members fReaderProxy or fValuePtr will be non-null
for a given TColumnValue, depending on whether the value comes from a real
TTree branch or from a temporary column respectively.

TDataFrame nodes can store tuples of TColumnValues and retrieve an updated
value for the column via the `Get` method.
**/
template <typename T>
class TColumnValue {
   // following line is equivalent to pseudo-code: ProxyParam_t == TArrayBranch<U> ? U : T
   // ReaderValueOrArray_t is a TTreeReaderValue<T> unless T is TArrayBranch<U>
   using ProxyParam_t = typename std::conditional<std::is_same<ReaderValueOrArray_t<T>, TTreeReaderValue<T>>::value, T,
                                                  TakeFirstParameter_t<T>>::type;

   /// TColumnValue has a slightly different behaviour whether the column comes from a TTreeReader, a TDataFrame Define
   /// or a TDataSource. It stores which it is as an enum.
   enum class EColumnKind { kTreeValue, kTreeArray, kCustomColumn, kDataSource, kInvalid };
   // Set to the correct value by MakeProxy or SetTmpColumn
   EColumnKind fColumnKind = EColumnKind::kInvalid;
   /// The slot this value belongs to. Only needed when querying custom column values, it is set in `SetTmpColumn`.
   unsigned int fSlot = std::numeric_limits<unsigned int>::max();

   // Each element of the following data members will be in use by a _single task_.
   // The vectors are used as very small stacks (1-2 elements typically) that fill in case of interleaved task execution
   // i.e. when more than one task needs readers in this worker thread.

   /// Owning ptrs to a TTreeReaderValue. Used for non-temporary columns when T != TArrayBranch<U>
   std::vector<std::unique_ptr<TTreeReaderValue<T>>> fReaderValues;
   /// Owning ptrs to a TTreeReaderArray. Used for non-temporary columns when T == TArrayBranch<U>.
   std::vector<std::unique_ptr<TTreeReaderArray<ProxyParam_t>>> fReaderArrays;
   /// Non-owning ptrs to the value of a custom column.
   std::vector<T *> fCustomValuePtrs;
   /// Non-owning ptrs to the value of a data-source column.
   std::vector<T **> fDSValuePtrs;
   /// Non-owning ptrs to the node responsible for the custom column. Needed when querying custom values.
   std::vector<TCustomColumnBase *> fCustomColumns;
   /// Signal whether we ever checked that the branch we are reading with a TTreeReaderArray stores array elements
   /// in contiguous memory. Only used when T == TArrayBranch<U>.
   bool fArrayHasBeenChecked = false;

public:
   TColumnValue() = default;

   void SetTmpColumn(unsigned int slot, TCustomColumnBase *tmpColumn);

   void MakeProxy(TTreeReader *r, const std::string &bn)
   {
      constexpr bool useReaderValue = std::is_same<ProxyParam_t, T>::value;
      if (useReaderValue) {
         fColumnKind = EColumnKind::kTreeValue;
         fReaderValues.emplace_back(new TTreeReaderValue<T>(*r, bn.c_str()));
      } else {
         fColumnKind = EColumnKind::kTreeArray;
         fReaderArrays.emplace_back(new TTreeReaderArray<ProxyParam_t>(*r, bn.c_str()));
      }
   }

   /// This overload is used to return scalar quantities (i.e. types that are not read into a TArrayBranch)
   template <typename U = T,
             typename std::enable_if<std::is_same<typename TColumnValue<U>::ProxyParam_t, U>::value, int>::type = 0>
   T &Get(Long64_t entry);

   /// This overload is used to return arrays (i.e. types that are read into a TArrayBranch)
   template <typename U = T, typename std::enable_if<!std::is_same<ProxyParam_t, U>::value, int>::type = 0>
   TArrayBranch<ProxyParam_t> Get(Long64_t)
   {
      auto &readerArray = *fReaderArrays.back();
      // We only use TTreeReaderArrays to read columns that users flagged as type `TArrayBranch`, so we need to check
      // that the branch stores the array as contiguous memory that we can actually wrap in an `TArrayBranch`.
      // Currently we need the first entry to have been loaded to perform the check
      // TODO Move check to `MakeProxy` once Axel implements this kind of check in TTreeReaderArray using TBranchProxy
      if (!fArrayHasBeenChecked) {
         if (readerArray.GetSize() > 1) {
            if (1 != (&readerArray[1] - &readerArray[0])) {
               std::string exceptionText = "Branch ";
               exceptionText += readerArray.GetBranchName();
               exceptionText +=
                  " hangs from a non-split branch. For this reason, it cannot be accessed via a TArrayBranch."
                  " Please read the top level branch instead.";
               throw std::runtime_error(exceptionText);
            }
            fArrayHasBeenChecked = true;
         }
      }

      return TArrayBranch<ProxyParam_t>(readerArray);
   }

   void Reset()
   {
      switch (fColumnKind) {
      case EColumnKind::kTreeValue: fReaderValues.pop_back(); break;
      case EColumnKind::kTreeArray: fReaderArrays.pop_back(); break;
      case EColumnKind::kCustomColumn:
         fCustomColumns.pop_back();
         fCustomValuePtrs.pop_back();
         break;
      case EColumnKind::kDataSource:
         fCustomColumns.pop_back();
         fDSValuePtrs.pop_back();
         break;
      case EColumnKind::kInvalid: throw std::runtime_error("ColumnKind not set for this TColumnValue");
      }
   }
};

template <typename T>
struct TTDFValueTuple {
};

template <typename... BranchTypes>
struct TTDFValueTuple<TypeList<BranchTypes...>> {
   using type = std::tuple<TColumnValue<BranchTypes>...>;
};

template <typename BranchType>
using TDFValueTuple_t = typename TTDFValueTuple<BranchType>::type;

/// Clear the proxies of a tuple of TColumnValues
template <typename ValueTuple, int... S>
void ResetTDFValueTuple(ValueTuple &values, StaticSeq<S...>)
{
   // hack to expand a parameter pack without c++17 fold expressions.
   std::initializer_list<int> expander{(std::get<S>(values).Reset(), 0)...};
   (void)expander; // avoid "unused variable" warnings
}

class TActionBase {
protected:
   TLoopManager *fImplPtr;     ///< A raw pointer to the TLoopManager at the root of this functional
                               /// graph. It is only guaranteed to contain a valid address during an
                               /// event loop.
   const unsigned int fNSlots; ///< Number of thread slots used by this node.

public:
   TActionBase(TLoopManager *implPtr, const unsigned int nSlots);
   TActionBase(const TActionBase &) = delete;
   TActionBase &operator=(const TActionBase &) = delete;
   virtual ~TActionBase() = default;

   virtual void Run(unsigned int slot, Long64_t entry) = 0;
   virtual void InitSlot(TTreeReader *r, unsigned int slot) = 0;
   virtual void TriggerChildrenCount() = 0;
   virtual void ClearValueReaders(unsigned int slot) = 0;
   unsigned int GetNSlots() const { return fNSlots; }
   /// This method is invoked to update a partial result during the event loop, right before passing the result to a
   /// user-defined callback registered via TResultProxy::RegisterCallback
   virtual void *PartialUpdate(unsigned int slot) = 0;
};

template <typename Helper, typename PrevDataFrame, typename BranchTypes_t = typename Helper::BranchTypes_t>
class TAction final : public TActionBase {
   using TypeInd_t = GenStaticSeq_t<BranchTypes_t::list_size>;

   Helper fHelper;
   const ColumnNames_t fBranches;
   PrevDataFrame &fPrevData;
   std::vector<TDFValueTuple_t<BranchTypes_t>> fValues;

public:
   TAction(Helper &&h, const ColumnNames_t &bl, PrevDataFrame &pd)
      : TActionBase(pd.GetImplPtr(), pd.GetNSlots()), fHelper(std::move(h)), fBranches(bl), fPrevData(pd),
        fValues(fNSlots)
   {
   }

   TAction(const TAction &) = delete;
   TAction &operator=(const TAction &) = delete;
   ~TAction() { fHelper.Finalize(); }

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      InitTDFValues(slot, fValues[slot], r, fBranches, fImplPtr->GetCustomColumnNames(), fImplPtr->GetBookedColumns(),
                    TypeInd_t());
      fHelper.InitSlot(r, slot);
   }

   void Run(unsigned int slot, Long64_t entry) final
   {
      // check if entry passes all filters
      if (fPrevData.CheckFilters(slot, entry))
         Exec(slot, entry, TypeInd_t());
   }

   template <int... S>
   void Exec(unsigned int slot, Long64_t entry, TDFInternal::StaticSeq<S...>)
   {
      (void)entry; // avoid bogus 'unused parameter' warning in gcc4.9
      fHelper.Exec(slot, std::get<S>(fValues[slot]).Get(entry)...);
   }

   void TriggerChildrenCount() final { fPrevData.IncrChildrenCount(); }

   virtual void ClearValueReaders(unsigned int slot) final { ResetTDFValueTuple(fValues[slot], TypeInd_t()); }

   /// This method is invoked to update a partial result during the event loop, right before passing the result to a
   /// user-defined callback registered via TResultProxy::RegisterCallback
   /// TODO the PartialUpdateImpl trick can go away once all action helpers will implement PartialUpdate
   void *PartialUpdate(unsigned int slot) final { return PartialUpdateImpl(slot); }

private:
   // this overload is SFINAE'd out if Helper does not implement `PartialUpdate`
   // the template parameter is required to defer instantiation of the method to SFINAE time
   template <typename H = Helper>
   auto PartialUpdateImpl(unsigned int slot) -> decltype(std::declval<H>().PartialUpdate(slot), (void *)(nullptr))
   {
      return &fHelper.PartialUpdate(slot);
   }
   // this one is always available but has lower precedence thanks to `...`
   void *PartialUpdateImpl(...) {
      throw std::runtime_error("This action does not support callbacks yet!");
   }
};

} // end NS TDF
} // end NS Internal

namespace Detail {
namespace TDF {

class TCustomColumnBase {
protected:
   TLoopManager *fImplPtr; ///< A raw pointer to the TLoopManager at the root of this functional graph. It is only
                           /// guaranteed to contain a valid address during an event loop.
   const std::string fName;
   unsigned int fNChildren{0};      ///< number of nodes of the functional graph hanging from this object
   unsigned int fNStopsReceived{0}; ///< number of times that a children node signaled to stop processing entries.
   const unsigned int fNSlots;      ///< number of thread slots used by this node, inherited from parent node.
   const bool fIsDataSourceColumn; ///< does the custom column refer to a data-source column? (or a user-define column?)
   std::vector<Long64_t> fLastCheckedEntry;

public:
   TCustomColumnBase(TLoopManager *df, std::string_view name, const unsigned int nSlots, const bool isDSColumn);
   TCustomColumnBase &operator=(const TCustomColumnBase &) = delete;
   virtual ~TCustomColumnBase() = default;

   virtual void InitSlot(TTreeReader *r, unsigned int slot) = 0;
   virtual void *GetValuePtr(unsigned int slot) = 0;
   virtual const std::type_info &GetTypeId() const = 0;
   TLoopManager *GetImplPtr() const;
   std::string GetName() const;
   virtual void Update(unsigned int slot, Long64_t entry) = 0;
   virtual void ClearValueReaders(unsigned int slot) = 0;
   unsigned int GetNSlots() const { return fNSlots; }
   bool IsDataSourceColumn() const { return fIsDataSourceColumn; }
   void InitNode();
};

namespace TCCHelperTypes {
struct TNothing {
};
struct TSlot {
};
struct TSlotAndEntry {
};
}

template <typename F, typename UPDATE_HELPER_TYPE = TCCHelperTypes::TNothing>
class TCustomColumn final : public TCustomColumnBase {
   // shortcuts
   using TNothing = TCCHelperTypes::TNothing;
   using TSlot = TCCHelperTypes::TSlot;
   using TSlotAndEntry = TCCHelperTypes::TSlotAndEntry;
   using UHT_t = UPDATE_HELPER_TYPE;
   // other types
   using FunParamTypes_t = typename CallableTraits<F>::arg_types;
   using BranchTypesTmp_t =
      typename TDFInternal::RemoveFirstParameterIf<std::is_same<TSlot, UHT_t>::value, FunParamTypes_t>::type;
   using BranchTypes_t = typename TDFInternal::RemoveFirstTwoParametersIf<std::is_same<TSlotAndEntry, UHT_t>::value,
                                                                          BranchTypesTmp_t>::type;
   using TypeInd_t = TDFInternal::GenStaticSeq_t<BranchTypes_t::list_size>;
   using ret_type = typename CallableTraits<F>::ret_type;
   // Avoid instantiating vector<bool> as `operator[]` returns temporaries in that case. Use std::deque instead.
   using ValuesPerSlot_t =
      typename std::conditional<std::is_same<ret_type, bool>::value, std::deque<ret_type>, std::vector<ret_type>>::type;

   F fExpression;
   const ColumnNames_t fBranches;
   ValuesPerSlot_t fLastResults;

   std::vector<TDFInternal::TDFValueTuple_t<BranchTypes_t>> fValues;

public:
   TCustomColumn(std::string_view name, F &&expression, const ColumnNames_t &bl, TLoopManager *lm,
                 bool isDSColumn = false)
      : TCustomColumnBase(lm, name, lm->GetNSlots(), isDSColumn), fExpression(std::move(expression)), fBranches(bl),
        fLastResults(fNSlots), fValues(fNSlots)
   {
   }

   TCustomColumn(const TCustomColumn &) = delete;
   TCustomColumn &operator=(const TCustomColumn &) = delete;

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      TDFInternal::InitTDFValues(slot, fValues[slot], r, fBranches, fImplPtr->GetCustomColumnNames(),
                                 fImplPtr->GetBookedColumns(), TypeInd_t());
   }

   void *GetValuePtr(unsigned int slot) final { return static_cast<void *>(&fLastResults[slot]); }

   void Update(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry[slot]) {
         // evaluate this filter, cache the result
         UpdateHelper(slot, entry, TypeInd_t(), BranchTypes_t(), (UPDATE_HELPER_TYPE *)nullptr);
         fLastCheckedEntry[slot] = entry;
      }
   }

   const std::type_info &GetTypeId() const
   {
      return fIsDataSourceColumn ? typeid(typename std::remove_pointer<ret_type>::type) : typeid(ret_type);
   }

   template <int... S, typename... BranchTypes>
   void UpdateHelper(unsigned int slot, Long64_t entry, TDFInternal::StaticSeq<S...>, TypeList<BranchTypes...>,
                     TCCHelperTypes::TNothing *)
   {
      fLastResults[slot] = fExpression(std::get<S>(fValues[slot]).Get(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   template <int... S, typename... BranchTypes>
   void UpdateHelper(unsigned int slot, Long64_t entry, TDFInternal::StaticSeq<S...>, TypeList<BranchTypes...>,
                     TCCHelperTypes::TSlot *)
   {
      fLastResults[slot] = fExpression(slot, std::get<S>(fValues[slot]).Get(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   template <int... S, typename... BranchTypes>
   void UpdateHelper(unsigned int slot, Long64_t entry, TDFInternal::StaticSeq<S...>, TypeList<BranchTypes...>,
                     TCCHelperTypes::TSlotAndEntry *)
   {
      fLastResults[slot] = fExpression(slot, entry, std::get<S>(fValues[slot]).Get(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   void ClearValueReaders(unsigned int slot) final { ResetTDFValueTuple(fValues[slot], TypeInd_t()); }
};

class TFilterBase {
protected:
   TLoopManager *fImplPtr; ///< A raw pointer to the TLoopManager at the root of this functional graph. It is only
                           /// guaranteed to contain a valid address during an event loop.
   std::vector<Long64_t> fLastCheckedEntry;
   std::vector<int> fLastResult = {true}; // std::vector<bool> cannot be used in a MT context safely
   std::vector<ULong64_t> fAccepted = {0};
   std::vector<ULong64_t> fRejected = {0};
   const std::string fName;
   unsigned int fNChildren{0};      ///< Number of nodes of the functional graph hanging from this object
   unsigned int fNStopsReceived{0}; ///< Number of times that a children node signaled to stop processing entries.
   const unsigned int fNSlots;      ///< Number of thread slots used by this node, inherited from parent node.

public:
   TFilterBase(TLoopManager *df, std::string_view name, const unsigned int nSlots);
   TFilterBase &operator=(const TFilterBase &) = delete;
   virtual ~TFilterBase() = default;

   virtual void InitSlot(TTreeReader *r, unsigned int slot) = 0;
   virtual bool CheckFilters(unsigned int slot, Long64_t entry) = 0;
   virtual void Report() const = 0;
   virtual void PartialReport() const = 0;
   TLoopManager *GetImplPtr() const;
   bool HasName() const;
   void PrintReport() const;
   virtual void IncrChildrenCount() = 0;
   virtual void StopProcessing() = 0;
   void ResetChildrenCount()
   {
      fNChildren = 0;
      fNStopsReceived = 0;
   }
   virtual void TriggerChildrenCount() = 0;
   unsigned int GetNSlots() const { return fNSlots; }
   void ResetReportCount()
   {
      assert(!fName.empty()); // this method is to only be called on named filters
      // fAccepted and fRejected could be different than 0 if this is not the first event-loop run using this filter
      std::fill(fAccepted.begin(), fAccepted.end(), 0);
      std::fill(fRejected.begin(), fRejected.end(), 0);
   }
   virtual void ClearValueReaders(unsigned int slot) = 0;
   void InitNode();
};

template <typename FilterF, typename PrevDataFrame>
class TFilter final : public TFilterBase {
   using BranchTypes_t = typename CallableTraits<FilterF>::arg_types;
   using TypeInd_t = TDFInternal::GenStaticSeq_t<BranchTypes_t::list_size>;

   FilterF fFilter;
   const ColumnNames_t fBranches;
   PrevDataFrame &fPrevData;
   std::vector<TDFInternal::TDFValueTuple_t<BranchTypes_t>> fValues;

public:
   TFilter(FilterF &&f, const ColumnNames_t &bl, PrevDataFrame &pd, std::string_view name = "")
      : TFilterBase(pd.GetImplPtr(), name, pd.GetNSlots()), fFilter(std::move(f)), fBranches(bl), fPrevData(pd),
        fValues(fNSlots)
   {
   }

   TFilter(const TFilter &) = delete;
   TFilter &operator=(const TFilter &) = delete;

   bool CheckFilters(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry[slot]) {
         if (!fPrevData.CheckFilters(slot, entry)) {
            // a filter upstream returned false, cache the result
            fLastResult[slot] = false;
         } else {
            // evaluate this filter, cache the result
            auto passed = CheckFilterHelper(slot, entry, TypeInd_t());
            passed ? ++fAccepted[slot] : ++fRejected[slot];
            fLastResult[slot] = passed;
         }
         fLastCheckedEntry[slot] = entry;
      }
      return fLastResult[slot];
   }

   template <int... S>
   bool CheckFilterHelper(unsigned int slot, Long64_t entry, TDFInternal::StaticSeq<S...>)
   {
      return fFilter(std::get<S>(fValues[slot]).Get(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      TDFInternal::InitTDFValues(slot, fValues[slot], r, fBranches, fImplPtr->GetCustomColumnNames(),
                                 fImplPtr->GetBookedColumns(), TypeInd_t());
   }

   // recursive chain of `Report`s
   void Report() const final { PartialReport(); }

   void PartialReport() const final
   {
      fPrevData.PartialReport();
      PrintReport();
   }

   void StopProcessing() final
   {
      ++fNStopsReceived;
      if (fNStopsReceived == fNChildren)
         fPrevData.StopProcessing();
   }

   void IncrChildrenCount() final
   {
      ++fNChildren;
      // propagate "children activation" upstream. named filters do the propagation via `TriggerChildrenCount`.
      if (fNChildren == 1 && fName.empty())
         fPrevData.IncrChildrenCount();
   }

   void TriggerChildrenCount() final
   {
      assert(!fName.empty()); // this method is to only be called on named filters
      fPrevData.IncrChildrenCount();
   }

   virtual void ClearValueReaders(unsigned int slot) final { ResetTDFValueTuple(fValues[slot], TypeInd_t()); }
};

class TRangeBase {
protected:
   TLoopManager *fImplPtr; ///< A raw pointer to the TLoopManager at the root of this functional graph. It is only
                           /// guaranteed to contain a valid address during an event loop.
   unsigned int fStart;
   unsigned int fStop;
   unsigned int fStride;
   Long64_t fLastCheckedEntry{-1};
   bool fLastResult{true};
   ULong64_t fNProcessedEntries{0};
   unsigned int fNChildren{0};      ///< Number of nodes of the functional graph hanging from this object
   unsigned int fNStopsReceived{0}; ///< Number of times that a children node signaled to stop processing entries.
   bool fHasStopped{false};         ///< True if the end of the range has been reached
   const unsigned int fNSlots;      ///< Number of thread slots used by this node, inherited from parent node.

public:
   TRangeBase(TLoopManager *implPtr, unsigned int start, unsigned int stop, unsigned int stride,
              const unsigned int nSlots);
   TRangeBase &operator=(const TRangeBase &) = delete;
   virtual ~TRangeBase() = default;

   TLoopManager *GetImplPtr() const;
   virtual bool CheckFilters(unsigned int slot, Long64_t entry) = 0;
   virtual void Report() const = 0;
   virtual void PartialReport() const = 0;
   virtual void IncrChildrenCount() = 0;
   virtual void StopProcessing() = 0;
   void ResetChildrenCount()
   {
      fNChildren = 0;
      fNStopsReceived = 0;
   }
   unsigned int GetNSlots() const { return fNSlots; }
};

template <typename PrevData>
class TRange final : public TRangeBase {
   PrevData &fPrevData;

public:
   TRange(unsigned int start, unsigned int stop, unsigned int stride, PrevData &pd)
      : TRangeBase(pd.GetImplPtr(), start, stop, stride, pd.GetNSlots()), fPrevData(pd)
   {
   }

   TRange(const TRange &) = delete;
   TRange &operator=(const TRange &) = delete;

   /// Ranges act as filters when it comes to selecting entries that downstream nodes should process
   bool CheckFilters(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry) {
         if (fHasStopped)
            return false;
         if (!fPrevData.CheckFilters(slot, entry)) {
            // a filter upstream returned false, cache the result
            fLastResult = false;
         } else {
            // apply range filter logic, cache the result
            ++fNProcessedEntries;
            if (fNProcessedEntries <= fStart || (fStop > 0 && fNProcessedEntries > fStop) ||
                (fStride != 1 && fNProcessedEntries % fStride != 0))
               fLastResult = false;
            else
               fLastResult = true;
            if (fNProcessedEntries == fStop) {
               fHasStopped = true;
               fPrevData.StopProcessing();
            }
         }
         fLastCheckedEntry = entry;
      }
      return fLastResult;
   }

   // recursive chain of `Report`s
   // TRange simply forwards these calls to the previous node
   void Report() const final { fPrevData.PartialReport(); }

   void PartialReport() const final { fPrevData.PartialReport(); }

   void StopProcessing() final
   {
      ++fNStopsReceived;
      if (fNStopsReceived == fNChildren && !fHasStopped)
         fPrevData.StopProcessing();
   }

   void IncrChildrenCount() final
   {
      ++fNChildren;
      // propagate "children activation" upstream
      if (fNChildren == 1)
         fPrevData.IncrChildrenCount();
   }
};

} // namespace TDF
} // namespace Detail

// method implementations
namespace Internal {
namespace TDF {

template <typename T>
void TColumnValue<T>::SetTmpColumn(unsigned int slot, ROOT::Detail::TDF::TCustomColumnBase *customColumn)
{
   fCustomColumns.emplace_back(customColumn);
   if (customColumn->GetTypeId() != typeid(T))
      throw std::runtime_error(
         std::string("TColumnValue: type specified for column \"" + customColumn->GetName() + "\" is ") +
         typeid(T).name() + " but temporary column has type " + customColumn->GetTypeId().name());

   if (customColumn->IsDataSourceColumn()) {
      fColumnKind = EColumnKind::kDataSource;
      fDSValuePtrs.emplace_back(static_cast<T **>(customColumn->GetValuePtr(slot)));
   } else {
      fColumnKind = EColumnKind::kCustomColumn;
      fCustomValuePtrs.emplace_back(static_cast<T *>(customColumn->GetValuePtr(slot)));
   }
   fSlot = slot;
}

// This method is executed inside the event-loop, many times per entry
// If need be, the if statement can be avoided using thunks
// (have both branches inside functions and have a pointer to
// the branch to be executed)
template <typename T>
template <typename U,
          typename std::enable_if<std::is_same<typename TColumnValue<U>::ProxyParam_t, U>::value, int>::type>
T &TColumnValue<T>::Get(Long64_t entry)
{
   if (fColumnKind == EColumnKind::kTreeValue) {
      return *(fReaderValues.back()->Get());
   } else {
      fCustomColumns.back()->Update(fSlot, entry);
      return fColumnKind == EColumnKind::kCustomColumn ? *fCustomValuePtrs.back() : **fDSValuePtrs.back();
   }
}

} // namespace TDF
} // namespace Internal
} // namespace ROOT
#endif // ROOT_TDFNODES
