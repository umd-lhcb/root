/*****************************************************************************

 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
// 
// BEGIN_HTML
// PiecewiseInterpolation 
// END_HTML
//

#include <exception>
#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooAbsReal.h"
#include "RooAbsSelfCachedReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooRealVar.h"
#include "RooMsgService.h"
#include "RooNumIntConfig.h"
#include "RooTrace.h"

using namespace std;

ClassImp(PiecewiseInterpolation)
;


//_____________________________________________________________________________
PiecewiseInterpolation::PiecewiseInterpolation()
{
  _positiveDefinite=false;
  TRACE_CREATE
}



//_____________________________________________________________________________
PiecewiseInterpolation::PiecewiseInterpolation(const char* name, const char* title, const RooAbsReal& nominal,
					       const RooArgList& lowSet, 
					       const RooArgList& highSet,
					       const RooArgList& paramSet,
					       Bool_t takeOwnership) :
  RooAbsSelfCachedReal(name, title, 0),
  _depList("!depList","list of observables", this),
  _nominal("!nominal","nominal value", this, (RooAbsReal&)nominal),
  _lowSet("!lowSet","low-side variation",this),
  _highSet("!highSet","high-side variation",this),
  _paramSet("!paramSet","high-side variation",this),
  _positiveDefinite(false)

{
  _depList.add(*(_nominal.arg().getObservables((((RooHistFunc&)_nominal.arg()).dataHist()))));
  // Constructor with two set of RooAbsReals. The value of the function will be
  //
  //  A = sum_i lowSet(i)*highSet(i) 
  //
  // If takeOwnership is true the PiecewiseInterpolation object will take ownership of the arguments in sumSet

  // KC: check both sizes
  if (lowSet.getSize() != highSet.getSize()) {
    coutE(InputArguments) << "PiecewiseInterpolation::ctor(" << GetName() << ") ERROR: input lists should be of equal length" << endl ;
    RooErrorHandler::softAbort() ;    
  }

  RooFIter inputIter1 = lowSet.fwdIterator() ;
  RooAbsArg* comp ;
  while((comp = inputIter1.next())) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "PiecewiseInterpolation::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _lowSet.add(*comp) ;
    if (takeOwnership) {
      _ownedList.addOwned(*comp) ;
    }
  }


  RooFIter inputIter2 = highSet.fwdIterator() ;
  while((comp = inputIter2.next())) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "PiecewiseInterpolation::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _highSet.add(*comp) ;
    if (takeOwnership) {
      _ownedList.addOwned(*comp) ;
    }
  }


  RooFIter inputIter3 = paramSet.fwdIterator() ;
  while((comp = inputIter3.next())) {
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      coutE(InputArguments) << "PiecewiseInterpolation::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " in first list is not of type RooAbsReal" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _paramSet.add(*comp) ;
    if (takeOwnership) {
      _ownedList.addOwned(*comp) ;
    }
    _interpCode.push_back(0); // default code: linear interpolation
  }

  
  // Choose special integrator by default 
  specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator") ;
  TRACE_CREATE
  this->Print("V");
}



//_____________________________________________________________________________
PiecewiseInterpolation::PiecewiseInterpolation(const PiecewiseInterpolation& other, const char* name) :
  RooAbsSelfCachedReal(other, name), 
  _depList("!depList",this,other._depList),
  _nominal("!nominal",this,other._nominal),
  _lowSet("!lowSet",this,other._lowSet),
  _highSet("!highSet",this,other._highSet),
  _paramSet("!paramSet",this,other._paramSet),
  _positiveDefinite(other._positiveDefinite),
  _interpCode(other._interpCode)
{
  // Copy constructor

  // Member _ownedList is intentionally not copy-constructed -- ownership is not transferred
  TRACE_CREATE
}



//_____________________________________________________________________________
PiecewiseInterpolation::~PiecewiseInterpolation() 
{
  // Destructor
  TRACE_DESTROY
}

RooArgSet* PiecewiseInterpolation::actualObservables(const RooArgSet &nset) const
{
  RooArgSet *myDeps = new RooArgSet;
  TIterator* diter = _depList.createIterator();
  diter->Reset();
  RooAbsArg* dep;
  while((dep=(RooAbsArg*)diter->Next()))
  {
    myDeps->add(*dep);
  }
  return myDeps;
  //return (RooArgSet*) myDeps.selectCommon(nset);


}

RooArgSet* PiecewiseInterpolation::actualParameters(const RooArgSet &nset) const
{
  RooArgSet *myPars=new RooArgSet;
  TIterator* piter = _paramSet.createIterator();
  piter->Reset();
  RooAbsArg* par;
  while((par=(RooAbsArg*)piter->Next()))
  {
    myPars->add(*par);
  }
  myPars->remove(nset,kTRUE,kTRUE);
  return myPars;


}


Double_t PiecewiseInterpolation::getValV(const RooArgSet* nset) const 
{
  return RooAbsCachedReal::getValV(((RooArgSet*)&_depList));


}
//8*/
void PiecewiseInterpolation::fillCacheObject(RooAbsCachedReal::FuncCacheElem& cache) const
{
  RooDataHist& cacheHist = *(cache.hist()) ;
  RooDataHist& nomHist=((RooHistFunc&)(_nominal.arg())).dataHist();


  // Iterator over all bins of RooDataHist and fill weights
  for (Int_t k=0 ; k<cacheHist.numEntries() ; k++) {
    cacheHist.get(k) ;

    nomHist.get(k);


    ///////////////////
    Double_t nominal = nomHist.weight();
    Double_t sum(nominal) ;

  //  return sum;
    RooAbsReal* param ;
    RooHistFunc* high ;
    RooHistFunc* low ;
    int i=-1;

    RooFIter lowIter(_lowSet.fwdIterator()) ;
    RooFIter highIter(_highSet.fwdIterator()) ;
    RooFIter paramIter(_paramSet.fwdIterator()) ;

    while((param=(RooAbsReal*)paramIter.next())) {
      ++i;
      low = (RooHistFunc*)(lowIter.next()) ;
      high = (RooHistFunc*)(highIter.next()) ;
      low->dataHist().get(k);
      high->dataHist().get(k);
      Double_t highVal=high->dataHist().weight();
      Double_t lowVal=low->dataHist().weight();
      if ((lowVal==nominal))
      {
        if(highVal==nominal) continue;
      }

      Int_t icode = _interpCode[i] ;

      switch(icode) {
      case 0: {
        // piece-wise linear
        if(param->getVal()>0)
    sum +=  param->getVal()*(highVal - nominal );
        else
    sum += param->getVal()*(nominal - lowVal);
        break ;
      }
      case 1: {
        // pice-wise log
        if(param->getVal()>=0)
    sum *= pow(highVal/nominal, +param->getVal());
        else
    sum *= pow(lowVal/nominal,  -param->getVal());
        break ;
      }
      case 2: {
        // parabolic with linear
        double a = 0.5*(highVal+lowVal)-nominal;
        double b = 0.5*(highVal-lowVal);
        double c = 0;
        if(param->getVal()>1 ){
    sum += (2*a+b)*(param->getVal()-1)+highVal-nominal;
        } else if(param->getVal()<-1 ) {
    sum += -1*(2*a-b)*(param->getVal()+1)+lowVal-nominal;
        } else {
    sum +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
        }
        break ;
      }
      case 3: {
        //parabolic version of log-normal
        double a = 0.5*(highVal+lowVal)-nominal;
        double b = 0.5*(highVal-lowVal);
        double c = 0;
        if(param->getVal()>1 ){
    sum += (2*a+b)*(param->getVal()-1)+highVal-nominal;
        } else if(param->getVal()<-1 ) {
    sum += -1*(2*a-b)*(param->getVal()+1)+lowVal-nominal;
        } else {
    sum +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
        }
        break ;
      }
      case 4: {
        
        // WVE ****************************************************************
        // WVE *** THIS CODE IS CRITICAL TO HISTFACTORY FIT CPU PERFORMANCE ***
        // WVE *** Do not modify unless you know what you are doing...      ***
        // WVE ****************************************************************

        double x  = param->getVal();      
        if (x>1) {
    sum += x*(highVal - nominal );
        } else if (x<-1) {
    sum += x*(nominal - lowVal);
        } else {
    double eps_plus = highVal - nominal;
    double eps_minus = nominal - lowVal;
    double S = 0.5 * (eps_plus + eps_minus);
    double A = 0.0625 * (eps_plus - eps_minus);
    
    //fcns+der+2nd_der are eq at bd

          double val = nominal + x * (S + x * A * ( 15 + x * x * (-10 + x * x * 3  ) ) ); 


    if (val < 0) val = 0;
    sum += val-nominal;
        }
        break ;

        // WVE ****************************************************************
      }
      case 5: {
        
        double x0 = 1.0;//boundary;
        double x  = param->getVal();

        if (x > x0 || x < -x0)
        {
    if(x>0)
      sum += x*(highVal - nominal );
    else
      sum += x*(nominal - lowVal);
        }
        else if (nominal != 0)
        {
    double eps_plus = highVal - nominal;
    double eps_minus = nominal - lowVal;
    double S = (eps_plus + eps_minus)/2;
    double A = (eps_plus - eps_minus)/2;

    //fcns+der are eq at bd
    double a = S;
    double b = 3*A/(2*x0);
    //double c = 0;
    double d = -A/(2*x0*x0*x0);

    double val = nominal + a*x + b*pow(x, 2) + 0/*c*pow(x, 3)*/ + d*pow(x, 4);
    if (val < 0) val = 0;

    //cout << "Using interp code 5, val = " << val << endl;

    sum += val-nominal;
        }
        break ;
      }
      default: {
        coutE(InputArguments) << "PiecewiseInterpolation::evaluate ERROR:  " << param->GetName() 
            << " with unknown interpolation code" << icode << endl ;
        break ;
      }
      }
    }
    
    if(_positiveDefinite && (sum<0)){
      sum = 1e-6;
      sum = 0;
      //     cout <<"sum < 0 forcing  positive definite"<<endl;
       //     int code = 1;
       //     RooArgSet* myset = new RooArgSet();
       //     cout << "integral = " << analyticalIntegralWN(code, myset) << endl;
    } else if(sum<0){
       cxcoutD(Tracing) <<"PiecewiseInterpolation::evaluate -  sum < 0, not forcing positive definite"<<endl;
    }
    cacheHist.set(sum);

  }

}


//_____________________________________________________________________________
Double_t PiecewiseInterpolation::evaluate() const 
{
  // Calculate and return current value of self

  ///////////////////
  Double_t nominal = _nominal;
  Double_t sum(nominal) ;

//  return sum;
  RooAbsReal* param ;
  RooAbsReal* high ;
  RooAbsReal* low ;
  int i=-1;

  RooFIter lowIter(_lowSet.fwdIterator()) ;
  RooFIter highIter(_highSet.fwdIterator()) ;
  RooFIter paramIter(_paramSet.fwdIterator()) ;

  while((param=(RooAbsReal*)paramIter.next())) {
    ++i;
    low = (RooAbsReal*)lowIter.next() ;
    high = (RooAbsReal*)highIter.next() ;
    Double_t highVal=high->getVal();
    Double_t lowVal=low->getVal();
    if ((lowVal==nominal))
    {
      if(highVal==nominal) continue;
    }

    Int_t icode = _interpCode[i] ;

    switch(icode) {
    case 0: {
      // piece-wise linear
      if(param->getVal()>0)
	sum +=  param->getVal()*(highVal - nominal );
      else
	sum += param->getVal()*(nominal - lowVal);
      break ;
    }
    case 1: {
      // pice-wise log
      if(param->getVal()>=0)
	sum *= pow(highVal/nominal, +param->getVal());
      else
	sum *= pow(lowVal/nominal,  -param->getVal());
      break ;
    }
    case 2: {
      // parabolic with linear
      double a = 0.5*(highVal+lowVal)-nominal;
      double b = 0.5*(highVal-lowVal);
      double c = 0;
      if(param->getVal()>1 ){
	sum += (2*a+b)*(param->getVal()-1)+highVal-nominal;
      } else if(param->getVal()<-1 ) {
	sum += -1*(2*a-b)*(param->getVal()+1)+lowVal-nominal;
      } else {
	sum +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
      }
      break ;
    }
    case 3: {
      //parabolic version of log-normal
      double a = 0.5*(highVal+lowVal)-nominal;
      double b = 0.5*(highVal-lowVal);
      double c = 0;
      if(param->getVal()>1 ){
	sum += (2*a+b)*(param->getVal()-1)+highVal-nominal;
      } else if(param->getVal()<-1 ) {
	sum += -1*(2*a-b)*(param->getVal()+1)+lowVal-nominal;
      } else {
	sum +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
      }
      break ;
    }
    case 4: {
      
      // WVE ****************************************************************
      // WVE *** THIS CODE IS CRITICAL TO HISTFACTORY FIT CPU PERFORMANCE ***
      // WVE *** Do not modify unless you know what you are doing...      ***
      // WVE ****************************************************************

      double x  = param->getVal();      
      if (x>1) {
	sum += x*(highVal - nominal );
      } else if (x<-1) {
	sum += x*(nominal - lowVal);
      } else {
	double eps_plus = highVal - nominal;
	double eps_minus = nominal - lowVal;
	double S = 0.5 * (eps_plus + eps_minus);
	double A = 0.0625 * (eps_plus - eps_minus);
	
	//fcns+der+2nd_der are eq at bd

        double val = nominal + x * (S + x * A * ( 15 + x * x * (-10 + x * x * 3  ) ) ); 


	if (val < 0) val = 0;
	sum += val-nominal;
      }
      break ;

      // WVE ****************************************************************
    }
    case 5: {
      
      double x0 = 1.0;//boundary;
      double x  = param->getVal();

      if (x > x0 || x < -x0)
      {
	if(x>0)
	  sum += x*(highVal - nominal );
	else
	  sum += x*(nominal - lowVal);
      }
      else if (nominal != 0)
      {
	double eps_plus = highVal - nominal;
	double eps_minus = nominal - lowVal;
	double S = (eps_plus + eps_minus)/2;
	double A = (eps_plus - eps_minus)/2;

	//fcns+der are eq at bd
	double a = S;
	double b = 3*A/(2*x0);
	//double c = 0;
	double d = -A/(2*x0*x0*x0);

	double val = nominal + a*x + b*pow(x, 2) + 0/*c*pow(x, 3)*/ + d*pow(x, 4);
	if (val < 0) val = 0;

	//cout << "Using interp code 5, val = " << val << endl;

	sum += val-nominal;
      }
      break ;
    }
    default: {
      coutE(InputArguments) << "PiecewiseInterpolation::evaluate ERROR:  " << param->GetName() 
			    << " with unknown interpolation code" << icode << endl ;
      break ;
    }
    }
  }
  
  if(_positiveDefinite && (sum<0)){
    sum = 1e-6;
    sum = 0;
    //     cout <<"sum < 0 forcing  positive definite"<<endl;
     //     int code = 1;
     //     RooArgSet* myset = new RooArgSet();
     //     cout << "integral = " << analyticalIntegralWN(code, myset) << endl;
  } else if(sum<0){
     cxcoutD(Tracing) <<"PiecewiseInterpolation::evaluate -  sum < 0, not forcing positive definite"<<endl;
  }
  return sum;

}

//_____________________________________________________________________________
Bool_t PiecewiseInterpolation::setBinIntegrator(RooArgSet& allVars) 
{

  if(allVars.getSize()==1){
    RooAbsReal* temp = const_cast<PiecewiseInterpolation*>(this);
    temp->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator")  ;
    int nbins = ((RooRealVar*) allVars.first())->numBins();
    temp->specialIntegratorConfig(kTRUE)->getConfigSection("RooBinIntegrator").setRealValue("numBins",nbins);
    return true;
  }else{
    cout << "Currently BinIntegrator only knows how to deal with 1-d "<<endl;
    return false;
  }
  return false;
}

//_____________________________________________________________________________
  
/*
  // KC: check if interCode=0 for all 
  RooFIter paramIterExtra(_paramSet.fwdIterator()) ;
  int i=0;
  while( paramIterExtra.next() ) {
    if(!_interpCode.empty() && _interpCode[i]!=0){
      // can't factorize integral
      cout <<"can't factorize integral"<<endl;
      return 0;
    }
    ++i;
  }

  // Select subset of allVars that are actual dependents
  analVars.add(allVars) ;
  //  RooArgSet* normSet = normSet2 ? getObservables(normSet2) : 0 ;
  //  RooArgSet* normSet = getObservables();
  //  RooArgSet* normSet = 0;


  // Check if this configuration was created before
  Int_t sterileIdx(-1) ;
  CacheElem* cache = (CacheElem*) _normIntMgr.getObj(normSet,&analVars,&sterileIdx) ;
  if (cache) {
    return _normIntMgr.lastIndex()+1 ;
  }
  
  // Create new cache element
  cache = new CacheElem ;

  // Make list of function projection and normalization integrals 
  RooAbsReal *func ;
  //  const RooArgSet* nset = _paramList.nset() ;

  // do nominal
  func = (RooAbsReal*)(&_nominal.arg()) ;
  RooAbsReal* funcInt = func->createIntegral(analVars) ;
  cache->_funcIntList.addOwned(*funcInt) ;

  // do variations
  RooFIter lowIter(_lowSet.fwdIterator()) ;
  RooFIter highIter(_highSet.fwdIterator()) ;
  RooFIter paramIter(_paramSet.fwdIterator()) ;

  //  int i=0;
  i=0;
  while(paramIter.next() ) {
    func = (RooAbsReal*)lowIter.next() ;
    funcInt = func->createIntegral(analVars) ;
    cache->_lowIntList.addOwned(*funcInt) ;

    func = (RooAbsReal*)highIter.next() ;
    funcInt = func->createIntegral(analVars) ;
    cache->_highIntList.addOwned(*funcInt) ;
    ++i;
  }

  // Store cache element
  Int_t code = _normIntMgr.setObj(normSet,&analVars,(RooAbsCacheElement*)cache,0) ;

  return code+1 ; 
}




//_____________________________________________________________________________
Double_t PiecewiseInterpolation::analyticalIntegralWN(Int_t code, const RooArgSet* ,const char* ) const 
{
  // Implement analytical integrations by doing appropriate weighting from  component integrals
  // functions to integrators of components



  

  // old integral, only works for linear and not positive definite
  CacheElem* cache = (CacheElem*) _normIntMgr.getObjByIndex(code-1) ;
  if( cache==NULL ) {
    std::cout << "Error: Cache Element is NULL" << std::endl;
    throw std::exception();
  }

  // old integral, only works for linear and not positive definite
  RooFIter funcIntIter = cache->_funcIntList.fwdIterator() ;
  RooFIter lowIntIter = cache->_lowIntList.fwdIterator() ;
  RooFIter highIntIter = cache->_highIntList.fwdIterator() ;
  RooAbsReal *funcInt(0), *low(0), *high(0), *param(0) ;
  Double_t value(0) ;
  Double_t nominal(0);

  // get nominal 
  int i=0;
  while(( funcInt = (RooAbsReal*)funcIntIter.next())) {
    value += funcInt->getVal() ;
    nominal = value;
    i++;
  }
  if(i==0 || i>1) { cout << "problem, wrong number of nominal functions"<<endl; }

  // now get low/high variations
  i = 0;
  RooFIter paramIter(_paramSet.fwdIterator()) ;

  // KC: old interp code with new iterator
  while( (param=(RooAbsReal*)paramIter.next()) ) {
    low = (RooAbsReal*)lowIntIter.next() ;
    high = (RooAbsReal*)highIntIter.next() ;
    
    if(param->getVal()>0) {
      value += param->getVal()*(high->getVal() - nominal );
    } else {
      value += param->getVal()*(nominal - low->getVal());
    }
    ++i;
  }



  //  cout << "value = " << value <<endl;
  return value;
}
*/

//_____________________________________________________________________________
void PiecewiseInterpolation::setInterpCode(RooAbsReal& param, int code){

  int index = _paramSet.index(&param);
  if(index<0){
      coutE(InputArguments) << "PiecewiseInterpolation::setInterpCode ERROR:  " << param.GetName() 
			    << " is not in list" << endl ;
  } else {
      coutW(InputArguments) << "PiecewiseInterpolation::setInterpCode :  " << param.GetName() 
			    << " is now " << code << endl ;
    _interpCode.at(index) = code;
  }
}


//_____________________________________________________________________________
void PiecewiseInterpolation::setAllInterpCodes(int code){

  for(unsigned int i=0; i<_interpCode.size(); ++i){
    _interpCode.at(i) = code;
  }
}


//_____________________________________________________________________________
void PiecewiseInterpolation::printAllInterpCodes(){

  for(unsigned int i=0; i<_interpCode.size(); ++i){
    coutI(InputArguments) <<"interp code for " << _paramSet.at(i)->GetName() << " = " << _interpCode.at(i) <<endl;
  }
}


//_____________________________________________________________________________
std::list<Double_t>* PiecewiseInterpolation::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  // WVE note: assumes nominal and alternates have identical structure, must add explicit check
  return _nominal.arg().binBoundaries(obs,xlo,xhi) ;  
}


//_____________________________________________________________________________
Bool_t PiecewiseInterpolation::isBinnedDistribution(const RooArgSet& obs) const 
{
  // WVE note: assumes nominal and alternates have identical structure, must add explicit check
  return _nominal.arg().isBinnedDistribution(obs) ;
}



//_____________________________________________________________________________
std::list<Double_t>* PiecewiseInterpolation::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  return _nominal.arg().plotSamplingHint(obs,xlo,xhi) ;  
}

//______________________________________________________________________________
void PiecewiseInterpolation::Streamer(TBuffer &R__b)
{
   // Stream an object of class PiecewiseInterpolation.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PiecewiseInterpolation::Class(),this);
      specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator") ;      
      if (_interpCode.empty()) _interpCode.resize(_paramSet.getSize());
   } else {
      R__b.WriteClassBuffer(PiecewiseInterpolation::Class(),this);
   }
}


/*
//_____________________________________________________________________________
void PiecewiseInterpolation::printMetaArgs(ostream& os) const 
{
  // Customized printing of arguments of a PiecewiseInterpolation to more intuitively reflect the contents of the
  // product operator construction

  _lowIter->Reset() ;
  if (_highIter) {
    _highIter->Reset() ;
  }

  Bool_t first(kTRUE) ;
    
  RooAbsArg* arg1, *arg2 ;
  if (_highSet.getSize()!=0) { 

    while((arg1=(RooAbsArg*)_lowIter->Next())) {
      if (!first) {
	os << " + " ;
      } else {
	first = kFALSE ;
      }
      arg2=(RooAbsArg*)_highIter->Next() ;
      os << arg1->GetName() << " * " << arg2->GetName() ;
    }

  } else {
    
    while((arg1=(RooAbsArg*)_lowIter->Next())) {
      if (!first) {
	os << " + " ;
      } else {
	first = kFALSE ;
      }
      os << arg1->GetName() ; 
    }  

  }

  os << " " ;    
}

*/
