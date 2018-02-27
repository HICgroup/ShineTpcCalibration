#include "TMath.h"
#include "Math/SMatrix.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"
#include "TF1.h"
#include "TApplication.h"
#include <iostream>
#include <vector>
#include <map>
#include "TGraphAsymmErrors.h"
#include "TGraph.h"

#include <fstream>
#include <sstream>

#include "TGraphErrors.h"
#include "Math/ParamFunctor.h"



double startRange;
double stopRange;
const double maxInvalidEntryRate = 0.5;
const unsigned int minNEntries = 4;


enum eCommonPlanes{
  eMTPCLvsTOFL = 0,
  eMTPCRvsTOFR,
  eVTPC2vsMTPCL,
  eVTPC2vsMTPCR,
  eVTPC1vsVTPC2,
  eVTPC1vsGTPC,
  eGTPCvsVTPC2,
  eMTPCRvsVTPC2_crosscheck
};


double sqr(double value){
  return value*value;
}

class uncertain_double{
public:
  uncertain_double(){
    reset();
  }

  uncertain_double const & operator+=(double rhs)
  {
    nentries++;
    sum+=rhs;
    sqr_sum+=rhs*rhs;
    return *this;
  }
  double estimate()
  {
    return sum/(double)nentries;
  }
  double error()
  {
    return (sqr_sum - sum*sum/(double)nentries)/(double)(nentries*(nentries - 1));
  }
  void reset()
  {
    sum = sqr_sum = 0;
    nentries = 0;
  }

private:
  double sum;
  double sqr_sum;
  unsigned int nentries;
};

struct line{
  ROOT::Math::SVector<double, 2 > f;
  ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double,2> > Cov;

  std::vector<unsigned int> points;
  double chi2,sumw;
  int vdriftStartIndex; //the same is used for y0 shifts, dont be confused
  int vdriftStopIndex;
  //  chi2 = sum_i (y-ax-b)^2/sigmay^2 ;
  //  sum_i    w_i * (y^2   - 2 y*(ax + b) + (ax+b)^2 )
  //  0   =  sum_i  (-2yx + (ax+b)x)w
  //  0   =  sum_i  (-2y  + (ax+b))w
};


class fitter{
public:
  fitter(std::vector<double> const &ix,
	 std::vector<double> const &iy,
	 std::vector<double> const &iyError):
    x(ix),
    y(iy),
    yError(iyError)
  {
    sortedIndex.resize(x.size());
    TMath::Sort((unsigned int)x.size(), &x[0], &sortedIndex[0], false);
  }


  std::vector<line> const minimalExpectedChi2Fit(double const sigmaFactor){
    if(x.size() < 2) return std::vector<line>(0);
    std::vector<line> lines(1);
    lines[0].points.push_back( sortedIndex[0] );
    lines[0].points.push_back( sortedIndex[1] );
    fit( lines[0] );
   unsigned int lastUnusedIndex = 2;
   while(lastUnusedIndex<x.size()){
     std::cout << " minimalExpectedChi2Fit :  lastUnusedIndex = " << lastUnusedIndex << std::endl;
     line test = lines.back();
     test.points.push_back(sortedIndex[lastUnusedIndex]);
     fit(test);
     std::cout << " minimalExpectedChi2Fit : fit finished " << std::endl;
     //   if(test.chi2 + 2.0/test.sumw < lines.back().chi2 + 2.0/lines.back().sumw + 1.0/sqr(yError[sortedIndex[lastUnusedIndex]]))
     if(test.chi2 + 2.0*sigmaFactor  < lines.back().chi2 + 3.0*sigmaFactor)
       {
	 lines.back() = test;
	 lastUnusedIndex++;

	 std::cout << " minimalExpectedChi2Fit : point added to lines.back() " << std::endl;
       }
     else
       {
	 std::cout << " minimalExpectedChi2Fit : point used for new line " << std::endl;
	 if(lastUnusedIndex + 1 < x.size())
	   {
	     lines.push_back(line());
	     lines.back().points.push_back(sortedIndex[lastUnusedIndex++]);
	     lines.back().points.push_back(sortedIndex[lastUnusedIndex++]);
	     std::cout << " minimalExpectedChi2Fit : creating new line with two points, lastUnusedIndex =  " << lastUnusedIndex << std::endl;
	     fit(lines.back());
	     std::cout << " minimalExpectedChi2Fit : fit finished " << std::endl;

	   }
	 else{
	   std::cout << " minimalExpectedChi2Fit : last point used for constant line " << std::endl;
	   lines.push_back(line());
	   lines.back().points.push_back(sortedIndex[lastUnusedIndex]);
	   lines.back().f[0] = y[sortedIndex[lastUnusedIndex]];
	   lines.back().f[1] = 0;
	   lines.back().Cov(0,0) = sqr(yError[sortedIndex[lastUnusedIndex++]]);
	   lines.back().Cov(0,1) = 0;
	   lines.back().Cov(1,1) = 0;
	   std::cout << " minimalExpectedChi2Fit : breaking loop " << std::endl;
	   break;
 	 }
       }
   }
   std::cout << " minimalExpectedChi2Fit: lines.size() = " << lines.size()  << std::endl;
   return lines;
  }

  std::vector<line> const orderedFit(double maximalSignificance = 3.0)
  {
    if(x.size() < 2) return std::vector<line>(0);
    std::vector<line> lines(1);
    lines[0].points.push_back( sortedIndex[0] );
    lines[0].points.push_back( sortedIndex[1] );
    fit(lines[0]);
    unsigned int lastUnusedIndex = 2;

    while(lastUnusedIndex<x.size()){
      std::cout << "orderedFit : lastUnusedIndex = " << lastUnusedIndex << std::endl;
      double ySignificance = significance(lines.back(), sortedIndex[lastUnusedIndex]);
      //      std::cout << " ySignificance = " << ySignificance << std::endl;
      if(ySignificance < maximalSignificance )
	{
	  std::cout << "orderedFit : adding point to last line" << std::endl;
	  lines.back().points.push_back(sortedIndex[lastUnusedIndex]);
	  fit(lines.back());
	  lastUnusedIndex++;
	}
      else
	{
	  std::cout << "orderedFit : creating new line" << std::endl;
	  if(lastUnusedIndex + 1 < x.size())
	    {
	      std::cout << "orderedFit : creating new line from two points" << std::endl;
	      lines.push_back(line());
	      lines.back().points.push_back(sortedIndex[lastUnusedIndex++]);
	      lines.back().points.push_back(sortedIndex[lastUnusedIndex++]);
	      fit(lines.back());
	    }
	  else{
	    std::cout << "orderedFit : creating constant line from last point" << std::endl;
	    lines.push_back(line());
	    lines.back().points.push_back(sortedIndex[lastUnusedIndex]);
	    lines.back().f[0] = y[sortedIndex[lastUnusedIndex]];
	    lines.back().f[1] = 0;
	    lines.back().Cov(0,0) = sqr(yError[sortedIndex[lastUnusedIndex++]]);
	    lines.back().Cov(0,1) = 0;
	    lines.back().Cov(1,1) = 0;
	    break;}
	}
    }

    std::cout << " orderedFit: lines.size() = " << lines.size()  << std::endl;
    return lines;
  }

  double significance(line const &dataGroup, unsigned int index)
  {
    double const yExtrapolated = dataGroup.f[0] + dataGroup.f[1]*x[index];
    double const yFitError = dataGroup.Cov(0,0) + 2*dataGroup.Cov(0,1)*x[index] + dataGroup.Cov(1,1)*sqr(x[index]);
    //    std::cout << " x = " << x[index] << " yExtrapolated = " << yExtrapolated << " +- " << sqrt(yFitError) << " y = " << y[index] << " +- " << yError[index] << std::endl;
    return fabs(yExtrapolated - y[index])/sqrt(yFitError + sqr(yError[index]));
  }

  void fit(line &dataGroup)
  {
    ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double,2> > G;
    ROOT::Math::SVector<double, 2> h;

    for(int i=0; i!= dataGroup.points.size(); i++)
      {
	unsigned int const index = dataGroup.points[i];
	double const weight = 1.0/sqr(yError[index]);
	G(0,0) += weight;
	G(0,1) += weight * x[index];
	G(1,1) += weight * sqr(x[index]);

	h[0] += weight * y[index];
	h[1] += weight * x[index] * y[index];
      }

    double sum_weight = G(0,0);
    G *= 1.0/sum_weight;
    h *= 1.0/sum_weight;

    G.InvertChol();
    dataGroup.f = G*h;
    dataGroup.Cov = G;
    dataGroup.Cov *= 1.0/sum_weight;

    double chi2 = 0;
    for(int i=0; i!=  dataGroup.points.size(); i++)
      {
	unsigned int const index = dataGroup.points[i];
	double const weight = 1.0/sqr(yError[index]);
	chi2 += weight*sqr(y[index] - dataGroup.f[0] - dataGroup.f[1]*x[index] );
      }
    dataGroup.chi2 = chi2;
    dataGroup.sumw = sum_weight;
  }

  inline double sqr(double value){
    return value*value;
  }

private:
  std::vector<double> const &x;
  std::vector<double> const &y;
  std::vector<double> const &yError;
  std::vector<unsigned int> sortedIndex;
};


class lineWrapper{
public:
  lineWrapper(line & myLine) : fLine(myLine)
  {}
  ~lineWrapper(){};
  double eval(double *x, double *p){
    return fLine.f[0] + fLine.f[1]*x[0];
  }
private:
  line const fLine;
};

class multiLineWrapper{
public:
  multiLineWrapper(std::vector<line> const & lines, std::vector<double> const &x) : fLines(lines),fX(x) {}
  ~multiLineWrapper(){}
  double eval(double *x, double *p) {
    double yExtrapolation, yExtrapolationError;
    if(x[0] < fX[fLines.front().points.front()] ) {
      evalSingleLine(0, x[0], yExtrapolation, yExtrapolationError  );

      return yExtrapolation;
    }
    if(x[0] > fX[fLines.back().points.back()] ) {
      evalSingleLine(fLines.size()-1, x[0], yExtrapolation, yExtrapolationError  );

      return yExtrapolation;
    }
    unsigned int lineIndex = 0;
    bool interPolation = false;
    for(unsigned int i=0; i != fLines.size(); i++)
      {
	if( x[0] >= fX[fLines[i].points.front()] && x[0] <= fX[fLines[i].points.back()]  )
	  {
	    lineIndex = i;
	    interPolation = false;
	    break;
	  }
	if(x[0] > fX[fLines[i].points.back()] && i!=fLines.size()-1 && x[0] < fX[fLines[i+1].points.front()] )
	  {
	    lineIndex = i;
	    interPolation = true;
	    break;
	  }
      }


    if(lineIndex == fLines.size()) {
      //error!
      lineIndex = fLines.size()-1; interPolation = false;
    }

    evalSingleLine(lineIndex, x[0], yExtrapolation, yExtrapolationError  );
    if(!interPolation) {return yExtrapolation;}

    bool smoothInterpolation = true;
    if(smoothInterpolation) {
      double yExtrapolation2, yExtrapolationError2;
      evalSingleLine(lineIndex+1, x[0], yExtrapolation2, yExtrapolationError2  );

      double const weightx = (x[0] - fX[fLines[lineIndex].points.back()]) / ( fX[fLines[lineIndex+1].points.front()] - fX[fLines[lineIndex].points.back()]  );

      double const weight1 = (1.0-weightx)/yExtrapolationError;
      double const weight2 = (weightx)/yExtrapolationError2;
      double const sumw = weight1+weight2;
      double const yInterpolation = (yExtrapolation*weight1 + yExtrapolation2*weight2 )/sumw;
      return yInterpolation;
    }
    else {
      double y1,y2, yError1, yError2;
      evalSingleLine(lineIndex, fX[fLines[lineIndex].points.back()], y1, yError1  );
      evalSingleLine(lineIndex+1, fX[fLines[lineIndex+1].points.front()], y2, yError2  );
      double const weightx = (x[0] - fX[fLines[lineIndex].points.back()]) / ( fX[fLines[lineIndex+1].points.front()] - fX[fLines[lineIndex].points.back()]  );

      return y1 + (y2-y1)*weightx;
    }
  }

private:
  void evalSingleLine(unsigned int const lineIndex, double const x, double &y, double &yError)
  {
    y = fLines[lineIndex].f[0] + fLines[lineIndex].f[1]*x;
    yError = fLines[lineIndex].Cov(0,0) + 2*fLines[lineIndex].Cov(0,1)*x + fLines[lineIndex].Cov(1,1)*x*x;

  }

  std::vector<line> const fLines;
  std::vector<double > const fX;
};



/*
class weighted_uncertain_double{
  weighted_uncertain_double(){
    reset();
  }

  weighted_uncertain_double const & Add(double value, double sigma2){
    sum_value+= value/sigma2;
    sum_weight += 1.0/sigma2;
  }

  double estimate()
  {
    return sum_value/sum_weights;
  }

  double error()
  {
    return 1/sum_weights;
  }

};
*/

struct alignmentparams {
  std::string s_slave;
  std::string s_master;
  std::string s_TPCvsTPC;

  TH1F * h_slave, *h_master, *h_residual;
  //  TH2F * h_residual_vs_masterY;
  TH2F * h_residual_vs_slave;

  TGraphAsymmErrors * g_slope;
  TGraphAsymmErrors * g_offset;

  struct alignment_section {
    double slope;
    double offset;

    double cov_slope;
    double cov_offset;
    double cov_slope_offset;

    uncertain_double residual;
    uncertain_double master;
    uncertain_double slave;
    unsigned int starttime, endtime;
    unsigned int n_entries;
    unsigned int n_dropped_entries;

    void Print(){
      std::cout << "starttime " << starttime << " endtime " << endtime << std::endl;
      std::cout << "timediff " << endtime - starttime<< std::endl;
      std::cout << "slope = " << slope << " +- " << cov_slope << " offset = " << offset << " +- " << cov_offset << std::endl;

      std::cout << "master mean = " << master.estimate() << " cov = "<< master.error() << std::endl;
      std::cout << "slave mean = " << slave.estimate() << " cov = " <<slave.error() << std::endl;
      std::cout << "residual mean = " << residual.estimate() << " cov = " << residual.error() << std::endl;
    }
  };
  std::vector<alignment_section> sections;

  std::vector<double> recChamberXCenters;
  std::vector<double> recChamberYCenters;
  std::vector<double> recChamberZCenters;
  std::vector<double> recChamberRotAngleXs;
  std::vector<double> recChamberRotAngleYs;
  std::vector<double> recChamberRotAngleZs;
  std::vector<double> runNumbers;
};

struct vdrift_corrections{
  std::string tpcName;
  TGraphAsymmErrors * g_vdrift_correction;
  struct vdrift_section{
    double vdrift_scale;
    double vdrift_scale_error;
    unsigned int starttime, endtime;
  };
  std::vector<vdrift_section> sections;
};

struct y0_corrections{
  std::string tpcName;
  TGraphAsymmErrors * g_y0_correction;

  std::vector<line> linearFits; //only for y0 shifts, relative vdrift plots are not smoothed
  double factor; unsigned int timeOffset; //plots use relative time while files need absolute time
  multiLineWrapper * smoothFit;
  ROOT::Math::ParamFunctor * smoothFunctor;
  TF1 * smoothFunc;

  struct y_shift_section{
    double y_shift;
    double y_shift_error;
    unsigned int starttime, endtime;
  };
  std::vector<y_shift_section> sections;

  struct RecGeomParameters {
    unsigned int runNumberMin, runNumberMax;
    double chamberXCenter;
    double chamberYCenter;
    double chamberZCenter;
    double chamberRotAngleX;
    double chamberRotAngleY;
    double chamberRotAngleZ;
  };
  RecGeomParameters recGeomParameters;
};


class straightFitter{
public:
  straightFitter(){
    Reset();
  }
  void AddEntry(double x, double y)
  {
    meany+= y;
    meanx+= x;
    meanxy+=x*y;
    meanxx+=x*x;
    meanxyy+=x*y*y;
    nentries++;

    meanyy+=y*y;
    meanxxy+=x*x*y;
    meanxxyy+=x*x*y*y;
    meanxxx+=x*x*x;
    meanxxxx+=x*x*x*x;
    meanxxxy+=x*x*x*y;

  }

  void SetBias(double const offset,
	       double const slope,
	       double const offset_std,
	       double const slope_std,
	       double const offset_slope_cov )
  {
    bias_offset           = offset;
    bias_slope            = slope;
    bias_offset_std       = offset_std;
    bias_slope_std        = slope_std;
    bias_offset_slope_cov = offset_slope_cov;
  }

  bool GetResults(double &offset, double &slope, double &offset_std, double &slope_std, double &offset_slope_cov)
  {
    if(nentries<10) return false;
    double temp[4] = {1, meanx.estimate(), meanx.estimate(), meanxx.estimate()};
    //    ROOT::Math::SMatrix<double, 2> Ginv(temp, temp+4);
    ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double,2> > Ginv;
    Ginv(0,0) = temp[0];
    Ginv(0,1) = temp[1];
    Ginv(1,1) = temp[3];

    if( !Ginv.InvertChol() ) {std::cout << "ERROR : matrix was not positive definit" << std::endl; return false;}
    ROOT::Math::SVector<double, 2> h(meany.estimate(), meanxy.estimate());
    ROOT::Math::SVector<double, 2> F = Ginv*h;

    offset = F[0];
    slope =  F[1];

    double corr_meany_meanxy = (meanxyy.estimate() - meany.estimate()*meanxy.estimate())/(double)(nentries-1);
    double temp2[4] = { meany.error(), corr_meany_meanxy, corr_meany_meanxy, meanxy.error()  };
    ROOT::Math::SMatrix<double, 2> Cov_meany_meanxy(temp2, temp2+4);
    ROOT::Math::SMatrix<double, 2> Cov_offset_slope(Ginv * Cov_meany_meanxy * Ginv);

    std::cout << "Ginv : " << std::endl << Ginv << std::endl;
    std::cout << "Cov_meany_meanxy : " << std::endl << Cov_meany_meanxy << std::endl;
    std::cout << "Cov_offset_slope : " << std::endl << Cov_offset_slope << std::endl;


    double temp3[16] = { meanyy.estimate(),  meanxyy.estimate(),  meanxy.estimate(),  meanxxy.estimate(),
                         meanxyy.estimate(), meanxxyy.estimate(), meanxxy.estimate(), meanxxxy.estimate(),
                         meanxy.estimate(),  meanxxy.estimate(),  meanxx.estimate(),  meanxxx.estimate(),
                         meanxxy.estimate(), meanxxxy.estimate(), meanxxx.estimate(), meanxxxx.estimate()};
    ROOT::Math::SMatrix<double, 4> Cov_all_params(temp3, temp3+16);
    double temp4[4]  = { meany.estimate(), meanxy.estimate(), meanx.estimate(), meanxx.estimate() };

    for(int i=0; i!=4; i++)
      for(int j=0; j!=4; j++)
        Cov_all_params(i, j) -= temp4[i]*temp4[j];

    Cov_all_params*=1.0/(double)(nentries-1);

    double temp5[4]  = {0, 1, 1, 0};
    ROOT::Math::SMatrix<double, 2> dG_over_dx(temp5, temp5+4);
    double temp6[4]  = {0,0,0,1};
    ROOT::Math::SMatrix<double, 2> dG_over_dxx(temp6, temp6+4);
    ROOT::Math::SMatrix<double, 2> dGinv_over_dx(   -Ginv*dG_over_dx*Ginv) ;
    ROOT::Math::SMatrix<double, 2> dGinv_over_dxx( -Ginv*dG_over_dxx*Ginv) ;


    ROOT::Math::SVector<double, 2> dF_over_dx = dGinv_over_dx*h;
    ROOT::Math::SVector<double, 2> dF_over_dxx = dGinv_over_dxx*h;

    std::cout << "h : " << std::endl << h << std::endl;
    std::cout << "F : " << std::endl << F << std::endl;
    std::cout << "dGinv_over_dxx : " << std::endl << dGinv_over_dxx << std::endl;


    ROOT::Math::SMatrix<double, 2, 4 > dF_over_dparam;
    ROOT::Math::SMatrix<double, 4, 2 > dF_over_dparam_T;
    for(int i=0; i!=2; i++) {
      for(int j=0; j!=4; j++) {
        if(j==0 || j==1) {
          dF_over_dparam[i][j] = Ginv[i][j];
          dF_over_dparam_T[j][i] = Ginv[i][j];
        }
        if(j==2){
          dF_over_dparam[i][j] = dF_over_dx[i];
          dF_over_dparam_T[j][i] = dF_over_dx[i];
        }
        if(j==3){
          dF_over_dparam[i][j] = dF_over_dxx[i];
          dF_over_dparam_T[j][i] = dF_over_dxx[i];
        }
      }
    }

    ROOT::Math::SMatrix<double, 2 > Cov_offset_slope_full(dF_over_dparam * Cov_all_params * dF_over_dparam_T   );


    std::cout << "dF_over_dparam : " << std::endl << dF_over_dparam << std::endl;
    std::cout << "dF_over_dparam_T : " << std::endl << dF_over_dparam_T << std::endl;
    std::cout << "Cov_all_params: " << std::endl << Cov_all_params << std::endl;
    std::cout << "Cov_offset_slope_full : " << std::endl << Cov_offset_slope_full << std::endl;

    /// we don't need bootstrap to evaluate the fitting uncertainty
    offset_std = Cov_offset_slope_full(0,0);
    slope_std  = Cov_offset_slope_full(1,1);
    offset_slope_cov = Cov_offset_slope_full(0,1);

    /*offset_std = Cov_offset_slope(0,0);
    slope_std  = Cov_offset_slope(1,1);
    offset_slope_cov = Cov_offset_slope(0,1);
    */

    /*
    std::cout << " bias slope            = " << bias_slope << std::endl
              << " bias offset           = " << bias_offset << std::endl
              << " bias slope_std        = " << bias_slope_std << std::endl
              << " bias offset_std       = " << bias_offset_std << std::endl
              << " bias offset_slope_cov = " << bias_offset_slope_cov << std::endl;

    slope += bias_slope;
    offset+= bias_offset;
    slope_std += bias_slope_std;
    offset_std += bias_offset_std;
    offset_slope_cov += bias_offset_slope_cov;*/
    std::cout << "slope = " << slope << " +- " << sqrt(slope_std) << std::endl;
    return true;
  }

  void Reset()
  {
    meanx.reset();
    meanxy.reset();
    meanxx.reset();
    meany.reset();
    meanxyy.reset();

    meanyy.reset();
    meanxxy.reset();
    meanxxyy.reset();
    meanxxx.reset();
    meanxxxx.reset();
    meanxxxy.reset();
    nentries = 0;

    bias_offset = bias_slope = bias_offset_std = bias_slope_std = bias_offset_slope_cov = 0;
  }

  uncertain_double meany,meanxy,meanx,meanxx,meanxyy;
  uncertain_double meanyy, meanxxy,
		   meanxxyy, meanxxx,
		   meanxxxx, meanxxxy;
  unsigned int nentries;

private:
  double bias_offset;
  double bias_slope;
  double bias_offset_std;
  double bias_slope_std;
  double bias_offset_slope_cov;
};

struct vdriftparams{
  std::string tpcName;
  TH1F * h_slave, *h_master;
  //  TH2F * h_residual_vs_masterY;
  TH2F * h_master_vs_slave;

  TGraphAsymmErrors * g_slope;
  TGraphAsymmErrors * g_offset;
  TGraph * g_recVdrift;
  unsigned int valid_size;

  std::vector<line> linearFits;
  double factor; unsigned int timeOffset; //plots use relative time while files need absolute time
  multiLineWrapper * smoothFit;
  ROOT::Math::ParamFunctor * smoothFunctor;
  TF1 * smoothFunc;

  struct vdrift_section{
    vdrift_section(){validity = false;}
    double vdrift;
    double vdrift_error;
    double offset, cov_offset, cov_vdrift_offset;

    uncertain_double master_y;
    uncertain_double slave_t;
    unsigned int starttime, endtime;
    unsigned int n_entries;
    unsigned int n_dropped_entries;
    bool validity;
  };
  std::vector<vdrift_section> sections;
  int GetVdriftSection(unsigned int time) const
  {
    int ret = -1;
    for(int i=0; i!= sections.size(); i++)
      {
	//std::cout << "sections[i].starttime " << sections[i].starttime << " sections[i].endtime " << sections[i].endtime << std::endl;
	if(sections[i].starttime <= time && sections[i].endtime>=time){
	  if(sections[i].validity) {ret = i;} else {ret = -1;}
	  break;
	}
      }
    return ret;
  }

  std::vector<double> recVdrifts;
  std::vector<double> eventUnixTime;
};

void fitAbsVdrift(TTree * dataTree,
		  vdriftparams const & mastervdrift,
		  vdriftparams & slavevdrift,
		  bool swapMasterSlaveBranch = false,
		  unsigned int min_master_ndf = 40,
		  unsigned int min_slave_ndf = 40)
{
  bool needhistos = true;
  unsigned int eventSecond;
  unsigned int eventNumber;
  unsigned int runNumber;
  int section=0;
  double master, slave;

  double meany=0,meanxy=0, meanx=0, meanxx=0;
  double slave_cov_x = 0;
  double slave_cov_y = 0;
  double master_cov_x = -1;
  double master_cov_y = -1;
  double slave_chi2 = 0;
  double slave_ndf = 0;
  double master_chi2 = -1;
  double master_ndf = -1;
  unsigned int nbad_tracks =0;
  double master_vdrift = -1;
  double slave_vdrift = 0 ;
  double master_chamberxcenter = -1;
  double master_chamberycenter = -1;
  double master_chamberzcenter = -1;
  double master_chamberrotanglex = -1;
  double master_chamberrotangley = -1;
  double master_chamberrotanglez = -1;
  double slave_chamberxcenter = 0 ;
  double slave_chamberycenter = 0 ;
  double slave_chamberzcenter = 0 ;
  double slave_chamberrotanglex = 0 ;
  double slave_chamberrotangley = 0 ;
  double slave_chamberrotanglez = 0 ;
  if(dataTree == 0)
  {
    std::cout << "Tree is missing, skipping." << std::endl;
     return;
  }

  uncertain_double master_value, slave_value, residual;

  std::string s_slave = "slave";
  std::string s_master = "master";

  if(swapMasterSlaveBranch) {
    s_slave = "master";
    s_master = "slave";
  }

  dataTree->SetBranchAddress((s_slave + "_Y").c_str(), &slave );
  dataTree->SetBranchAddress((s_slave+"_cov_X_2").c_str(),  &slave_cov_x);
  dataTree->SetBranchAddress((s_slave+"_cov_Y_2").c_str(),  &slave_cov_y);
  dataTree->SetBranchAddress((s_slave+"_chi2").c_str(),  &slave_chi2);
  dataTree->SetBranchAddress((s_slave+"_ndf").c_str(),  &slave_ndf);
  dataTree->SetBranchAddress((s_slave+"_recVDrift").c_str(),  &slave_vdrift);
  dataTree->SetBranchAddress((s_slave+"_recChamberXCenter").c_str(),  &slave_chamberxcenter);
  dataTree->SetBranchAddress((s_slave+"_recChamberYCenter").c_str(),  &slave_chamberycenter);
  dataTree->SetBranchAddress((s_slave+"_recChamberZCenter").c_str(),  &slave_chamberzcenter);
  dataTree->SetBranchAddress((s_slave+"_recChamberRotAngleX").c_str(),  &slave_chamberrotanglex);
  dataTree->SetBranchAddress((s_slave+"_recChamberRotAngleY").c_str(),  &slave_chamberrotangley);
  dataTree->SetBranchAddress((s_slave+"_recChamberRotAngleZ").c_str(),  &slave_chamberrotanglez);

  bool tofSource = false;
  if(mastervdrift.tpcName == "TOFL" || mastervdrift.tpcName == "TOFR")
    {
      tofSource = true;
      dataTree->SetBranchAddress("tofY", &master );
    }
  else
    {
      dataTree->SetBranchAddress((s_master + "_Y").c_str(), &master );
      dataTree->SetBranchAddress((s_master+"_cov_X_2").c_str(),  &master_cov_x);
      dataTree->SetBranchAddress((s_master+"_cov_Y_2").c_str(),  &master_cov_y);
      dataTree->SetBranchAddress((s_master+"_chi2").c_str(),  &master_chi2);
      dataTree->SetBranchAddress((s_master+"_ndf").c_str(),  &master_ndf);
      dataTree->SetBranchAddress((s_master+"_recVDrift").c_str(),  &master_vdrift);
      dataTree->SetBranchAddress((s_master+"_recChamberXCenter").c_str(),  &master_chamberxcenter);
      dataTree->SetBranchAddress((s_master+"_recChamberYCenter").c_str(),  &master_chamberycenter);
      dataTree->SetBranchAddress((s_master+"_recChamberZCenter").c_str(),  &master_chamberzcenter);
      dataTree->SetBranchAddress((s_master+"_recChamberRotAngleX").c_str(),  &master_chamberrotanglex);
      dataTree->SetBranchAddress((s_master+"_recChamberRotAngleY").c_str(),  &master_chamberrotangley);
      dataTree->SetBranchAddress((s_master+"_recChamberRotAngleZ").c_str(),  &master_chamberrotanglez);
    }

  dataTree->SetBranchAddress("eventUnixTime", &eventSecond);
  dataTree->SetBranchAddress("eventNumber", &eventNumber);
  dataTree->SetBranchAddress("runNumber", &runNumber);
  dataTree->SetBranchAddress("sectionId", &section);

  if(needhistos == true)
    {
      double min_resy = -30000;
      double max_resy = 30000;
      double min_y = 0;
      double max_y = 0;
      std::string masterTPC = mastervdrift.tpcName;
      std::string slaveTPC = slavevdrift.tpcName;
      slavevdrift.h_master = new TH1F(("h_master" + masterTPC).c_str(), "master", 500, min_y, max_y);
      slavevdrift.h_master->GetXaxis()->SetTitle((masterTPC+"_Y [cm]").c_str());
      slavevdrift.h_slave = new TH1F(("h_slave" + slaveTPC).c_str(), "slave", 100, min_y, max_y);
      slavevdrift.h_slave->GetXaxis()->SetTitle((slaveTPC+" drift time [ns])").c_str());
      slavevdrift.h_master_vs_slave = new TH2F(("h_master_vs_slave_" + masterTPC + "_" + slaveTPC).c_str(),
                                                "master vs slave", 100, min_y, max_y, 500, min_resy, max_resy);
      slavevdrift.h_master_vs_slave->GetXaxis()->SetTitle((slaveTPC+" drift time [ns])").c_str() );
      slavevdrift.h_master_vs_slave->GetYaxis()->SetTitle((masterTPC+"_Y [cm]").c_str());

    }


  long int allnentries = dataTree->GetEntries();
  long int nentries =0;

  int oldsectionnumber = 0;
  unsigned int starttime = 0, oldstarttime = 0, endtime = 0;

  std::vector<std::vector<straightFitter> > fitter;
  std::vector<std::vector<int > > masterVdrifts;

  int old_master_vdrift_section_id = -1;

  std::cout << "FitAbsVdrift starts analysis" << std::endl;

  for(int i=0; i!= allnentries; i++) {
    dataTree->GetEntry(i);
    if(eventSecond<1000000000) {
      std::cout << "Skipping event from the 1980s" << std::endl;
      continue;
    }
    double slave_t  = slave;
    double master_t = master;

    slavevdrift.recVdrifts.push_back(slave_vdrift);
    std::cout << "slave_vdrift from init tree: " << slave_vdrift << std::endl;
    slavevdrift.eventUnixTime.push_back(eventSecond);
    std::cout << "eventSecond from init tree: " << eventSecond << std::endl;

    slave_t/=slave_vdrift;

    int const master_vdrift_section_id = mastervdrift.GetVdriftSection(eventSecond);
    if(!tofSource && master_vdrift_section_id!=-1) {
      master_t/=master_vdrift;
            //std::cout << "slave_vdrift = " << slave_vdrift << " master_vdrift = " << master_vdrift << " Cmaster_vdrift = " << mastervdrift.sections[master_vdrift_section_id].vdrift << std::endl;
    }
    /*    if(!tofSource && master_vdrift_section_id==-1) {
      //      std::cout << "!tofSource && master_vdrift_section_id==-1" << std::endl;
      old_master_vdrift_section_id = -1;
      continue;
      }*/


    if(section == -1)
    {
      std::cout << "unprocessed data, sectionId = -1" << std::endl;
    }

    if(section==0 && oldsectionnumber==0 ) {
      //      oldsectionnumber = section;
      continue;
    }

    if(section!=0 && section!=oldsectionnumber) {
      slavevdrift.sections.push_back(vdriftparams::vdrift_section());
      slavevdrift.sections.back().starttime = eventSecond;
      fitter.push_back(std::vector<straightFitter>(1));
      //      std::cout << "fitter.size() " << fitter.size() << " fitter.back().size() " << fitter.back().size() << std::endl;
      old_master_vdrift_section_id = master_vdrift_section_id;

      if(!tofSource) {
	masterVdrifts.push_back(std::vector<int>(1,  master_vdrift_section_id) );
      }
      slavevdrift.sections.back().n_entries=0;
      slavevdrift.sections.back().n_dropped_entries=0;
    }
    if(section!=0) {
      /// the endtime variable updated in every entry for this section, until the last entry is filled, so at the end it will contain the ending time
      slavevdrift.sections.back().endtime = eventSecond;
      slavevdrift.sections.back().n_entries++;

      //should read the correlation term too
      if( TMath::Sqrt( slave_cov_x + slave_cov_y   ) > 2
     	  || slave_ndf < min_slave_ndf
        || (  !tofSource && (  master_ndf < min_master_ndf || TMath::Sqrt( master_cov_x + master_cov_y  ) > 2 ) ) )
      {
	/*        parameters.sections.back().n_dropped_entries++;
                  std::cout << " uniqueName " << uniqueName
                  << " slave_cov_x " << slave_cov_x
                  << " slave_cov_y " << slave_cov_y
                  << " slave_ndf " << slave_ndf << " " << min_slave_ndf
                  << " master_cov_x " << master_cov_x
                  << " master_cov_y " << master_cov_y
                  << " master_ndf " << master_ndf << " " << min_master_ndf << std::endl;
	*/
        nbad_tracks++;
        oldsectionnumber = section;
        continue;
      }

      if(!tofSource && master_vdrift_section_id == -1 ) {
	if(masterVdrifts.back().size() == 1){
	  slavevdrift.sections.back().starttime = eventSecond;
	}
	continue;
      }
      //      std::cout << "good track " << old_master_vdrift_section_id << " " << master_vdrift_section_id << std::endl;
      if(!tofSource && !mastervdrift.sections[master_vdrift_section_id].validity) {
	std::cout << "xxx invalid master vdrift!?" << std::endl;
      }

      if(!tofSource && old_master_vdrift_section_id!=master_vdrift_section_id)
	{

	  masterVdrifts.back().push_back(master_vdrift_section_id);
	  fitter.back().push_back(straightFitter());

	}
      fitter.back().back().AddEntry(slave_t, master_t);


      slavevdrift.sections.back().master_y+=master_t;
      slavevdrift.sections.back().slave_t+=slave_t;

    }

    if(needhistos){
      slavevdrift.h_slave->Fill(slave_t);
      slavevdrift.h_master->Fill(master_t);
      slavevdrift.h_master_vs_slave->Fill(slave_t, master_t);
    }

    oldsectionnumber = section;
    old_master_vdrift_section_id = master_vdrift_section_id;
  }

  unsigned int size = slavevdrift.sections.size();
  int final_size = 0;
  std::cout << slavevdrift.tpcName << std::endl;
  for(int i=0; i!= size; i++)
    {
      double tempm[4] = {0,0,0,0};
      double tempv[2] = {0,0};
      //ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double,2> > sum_error_inv(tempm, tempm+4);
      ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double,2> > sum_error_inv; //numerical errors are significant if the matrix is not explicitly symmetric
      sum_error_inv(0,0) = 0;
      sum_error_inv(0,1) = 0;
      sum_error_inv(1,1) = 0;

      //ROOT::Math::SMatrix<double, 2 > sum_error_inv(tempm, tempm+4);
      ROOT::Math::SVector<double, 2 > best_fit(tempv, tempv+2);
      ROOT::Math::SVector<double, 2 > best_fit_2(tempv, tempv+2);
      int n_good_fits = 0;



      if(!tofSource){
	if(masterVdrifts[i].size()!=fitter[i].size()) { std::cout << "invalid masterVdrifts or fitter size" << std::endl;}
	for(int j=0; j!= fitter[i].size(); j++ ) {
	  //	  ROOT::Math::SMatrix< double, 2 > error(tempm, tempm+4);
	  ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double,2> > error;
	  error(0,0) = 0;
	  error(0,1) = 0;
	  error(1,1) = 0;
	  ROOT::Math::SVector< double, 2 > fit_params(tempv, tempv+2);
	  double slope, offset, cov_slope, cov_offset, cov_slope_offset;
	  slope = offset = cov_slope = cov_offset = cov_slope_offset = 0;
	  if(!mastervdrift.sections[masterVdrifts[i][j]].validity) {
	    std::cout << "invalid master vdrift section " << i << " " << j << std::endl;
	    continue;
	  }
	  if(fitter[i][j].GetResults(offset,
				     slope,
				     cov_offset,
				     cov_slope,
				     cov_slope_offset)){

	    double vdrift = slope * mastervdrift.sections[masterVdrifts[i][j]].vdrift;
	    double yoffset = offset*mastervdrift.sections[masterVdrifts[i][j]].vdrift
	      +mastervdrift.sections[masterVdrifts[i][j]].offset;
	    double vdrift_error = sqr(slope) * mastervdrift.sections[masterVdrifts[i][j]].vdrift_error
	      + sqr(cov_slope) * mastervdrift.sections[masterVdrifts[i][j]].vdrift;

	    double yoffset_error = mastervdrift.sections[masterVdrifts[i][j]].cov_offset
	      + sqr(offset)*mastervdrift.sections[masterVdrifts[i][j]].vdrift_error
	      + 2*offset*mastervdrift.sections[masterVdrifts[i][j]].cov_vdrift_offset
	      + sqr(mastervdrift.sections[masterVdrifts[i][j]].vdrift)*cov_offset;

	    double cov_yoffset_vdrift = cov_slope_offset*sqr(mastervdrift.sections[masterVdrifts[i][j]].vdrift)
	      + slope*mastervdrift.sections[masterVdrifts[i][j]].cov_vdrift_offset
	      + slope*mastervdrift.sections[masterVdrifts[i][j]].vdrift*mastervdrift.sections[masterVdrifts[i][j]].vdrift_error;

	    error(0,0) = yoffset_error;
	    error(0,1) = cov_yoffset_vdrift;
	    error(1,0) = error(0,1);
	    error(1,1) = vdrift_error;

	    fit_params[0] = yoffset;
	    fit_params[1] = vdrift;

	    bool invstatus  =  error.InvertChol();
	    std::cout << " inversion status : " << invstatus << std::endl;
	    if(invstatus) {
	      n_good_fits++;
	      best_fit+= error*fit_params;
	      //best_fit = fit_params;
	      sum_error_inv+=error;
	      std::cout << " vdrift = " << vdrift << " vdrift_error = " << vdrift_error << " yoffset = " << yoffset << " yoffset_error = " << yoffset_error << " cov_yoffset_vdrift = " << cov_yoffset_vdrift << std::endl;
	      std::cout << " error invert : " << std::endl << error << std::endl;
	    }
	  }
	  else {
	    std::cout << "Subregion skipped due to bad statistics" << std::endl;
	  }
	}
	if(n_good_fits>0) {
	  //ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double,2> > save_sum_error_inv(sum_error_inv);
	  //ROOT::Math::SMatrix< double, 2 > sum_error (sum_error_inv);
	  bool invstatus  = sum_error_inv.InvertChol();
	  std::cout << " inversion status : " << invstatus << std::endl;

	  //ROOT::Math::SMatrix< double, 2 > sum_error_invert(sum_error_inv);
	  //std::cout << " error*error_inv = " << std::endl << sum_error_invert*sum_error << std::endl << sum_error*sum_error_invert << std::endl;

	  std::cout << " best_fit : " << std::endl << best_fit << std::endl;
	  best_fit_2 = sum_error_inv*best_fit;
	  //best_fit_2 = best_fit;

	  slavevdrift.sections[i].offset = best_fit_2[0];
	  slavevdrift.sections[i].vdrift = best_fit_2[1];
	  slavevdrift.sections[i].cov_offset = sum_error_inv(0,0);
	  slavevdrift.sections[i].vdrift_error = sum_error_inv(1,1);
	  slavevdrift.sections[i].cov_vdrift_offset = sum_error_inv(1,0);
	  if(invstatus) {slavevdrift.sections[i].validity=true;}
	  else {slavevdrift.sections[i].validity=false;}


	  /*if(n_good_fits == 1) {slavevdrift.sections[i].validity=false;}
	  else {
	    std::cout << " n_good_fits = " << n_good_fits << " from " << fitter[i].size() << std::endl;
	    std::cout << " offset = " << best_fit_2[0] << " +- "  << sqrt(sum_error_inv(0,0)) << " sqrt(" << sum_error_inv(0,0) << ")" <<  std::endl;
	    std::cout << " slope = " << best_fit_2[1] << " +- "  << sqrt(sum_error_inv(1,1)) << " sqrt(" << sum_error_inv(1,1) << ")" << std::endl;
	    std::cout << " sum_error_inv : " << std::endl << sum_error_inv << std::endl;
	    if(slavevdrift.sections[i].vdrift < 0.00001 || slavevdrift.sections[i].vdrift > 1 ) {std::cout << "strange vdrift value " << std::endl;}
	    if(invstatus) final_size++;
	  }
	  */
	  if(invstatus) {final_size++;}


	} else {
	  slavevdrift.sections[i].validity=false;
	  std::cout << "No good subsection for this section" << std::endl;
	}
      } else {
		slavevdrift.sections[i].validity=fitter[i][0].GetResults(slavevdrift.sections[i].offset,
			      slavevdrift.sections[i].vdrift,
			      slavevdrift.sections[i].cov_offset,
			      slavevdrift.sections[i].vdrift_error,
			      slavevdrift.sections[i].cov_vdrift_offset);
		if(slavevdrift.sections[i].validity) {final_size++;}
      }



      /*for(int j=0; j!= fitter[i].size(); j++ )
	{
	  ROOT::Math::SMatrix<double, 2 > error(tempm, tempm+4);
	  ROOT::Math::SVector<double, 2 > fit_params(tempv, tempv+2);
	  if(fitter[i][j].GetResults(fit_params[0],
				     fit_params[1],
				     error(0,0),
				     error(1,1),
				     error(1,0)))
	    {
	      n_good_fits++;
	      error(0,1) = error(1,0);
	      error.Invert();
	      best_fit+= error*fit_params;
	      sum_error_inv+=error;
	    }
	  else
	    {
	      std::cout << "Subregion skipped due to bad statistics" << std::endl;
	    }
	}
      if(n_good_fits>0){
	sum_error_inv.Invert();
	best_fit = sum_error_inv*best_fit;
	slavevdrift.sections[i].offset = best_fit[0];
	slavevdrift.sections[i].vdrift = best_fit[1];
	slavevdrift.sections[i].cov_offset = sum_error_inv(0,0);
	slavevdrift.sections[i].vdrift_error = sum_error_inv(1,1);
	slavevdrift.sections[i].cov_vdrift_offset = sum_error_inv(1,0);

      }
      else{
	std::cout << "No good subsection for this section" << std::endl;
	}*/


      /*
      double slope, offset, cov_slope, cov_offset, cov_slope_offset;
      slope = offset = cov_slope = cov_offset = cov_slope_offset = 0;

      if(!tofSource){
	slavevdrift.sections[i].validity=
	  fitter[i][0].GetResults(offset,
				  slope,
				  cov_offset,
				  cov_slope,
				  cov_slope_offset);
	double vdrift = slope * mastervdrift.sections[masterVdrifts[i][0]].vdrift;
	double yoffset = offset*mastervdrift.sections[masterVdrifts[i][0]].vdrift + mastervdrift.sections[masterVdrifts[i][0]].offset;
	double vdrift_error = sqr(slope) * mastervdrift.sections[masterVdrifts[i][0]].vdrift_error +
	  sqr(cov_slope) * mastervdrift.sections[masterVdrifts[i][0]].vdrift;

	double yoffset_error = mastervdrift.sections[masterVdrifts[i][0]].cov_offset +
	  sqr(offset)*mastervdrift.sections[masterVdrifts[i][0]].vdrift_error+
	  2*offset*mastervdrift.sections[masterVdrifts[i][0]].cov_vdrift_offset+
	  sqr(mastervdrift.sections[masterVdrifts[i][0]].vdrift)*cov_offset;

	double cov_yoffset_vdrift = cov_slope_offset*sqr(mastervdrift.sections[masterVdrifts[i][0]].vdrift)
	  + slope*mastervdrift.sections[masterVdrifts[i][0]].cov_vdrift_offset
	  + slope*mastervdrift.sections[masterVdrifts[i][0]].vdrift*mastervdrift.sections[masterVdrifts[i][0]].vdrift_error;

	slavevdrift.sections[i].offset = yoffset;
	slavevdrift.sections[i].vdrift = vdrift;
	slavevdrift.sections[i].cov_offset = yoffset_error;
	slavevdrift.sections[i].vdrift_error = vdrift_error;
	slavevdrift.sections[i].cov_vdrift_offset = cov_yoffset_vdrift;
      }
      else {
	slavevdrift.sections[i].validity=fitter[i][0].GetResults(slavevdrift.sections[i].offset,
			      slavevdrift.sections[i].vdrift,
			      slavevdrift.sections[i].cov_offset,
			      slavevdrift.sections[i].vdrift_error,
			      slavevdrift.sections[i].cov_vdrift_offset);
      }
      */
    }
  slavevdrift.valid_size = final_size;
  std::cout << "FitAbsVdrift finished" << std::endl;
}

void fitVdrift(TTree * dataTree, alignmentparams & parameters ,
               bool needhistos,
               std::string uniqueName,
               bool useRecVdrift = false,
               unsigned int min_master_ndf = 40,
               unsigned int min_slave_ndf = 40)
{
  unsigned int eventSecond;
  unsigned int eventNumber;
  unsigned int runNumber;
  int section=0;
  double master, slave;
  //double dmasterY, dextrapY;
  double meany=0,meanxy=0, meanx=0, meanxx=0;
  double slave_cov_x = 0;
  double slave_cov_y = 0;
  double master_cov_x = -1;
  double master_cov_y = -1;
  double slave_chi2 = 0;
  double slave_ndf = 0;
  double master_chi2 = -1;
  double master_ndf = -1;
  unsigned int nbad_tracks =0;
  double master_vdrift = -1;
  double slave_vdrift = 0 ;
  double master_chamberxcenter = -1;
  double master_chamberycenter = -1;
  double master_chamberzcenter = -1;
  double master_chamberrotanglex = -1;
  double master_chamberrotangley = -1;
  double master_chamberrotanglez = -1;
  double slave_chamberxcenter = 0 ;
  double slave_chamberycenter = 0 ;
  double slave_chamberzcenter = 0 ;
  double slave_chamberrotanglex = 0 ;
  double slave_chamberrotangley = 0 ;
  double slave_chamberrotanglez = 0 ;

  parameters.s_TPCvsTPC = uniqueName;

  std::cout << "uniqueName " << uniqueName <<  " dataTree" << dataTree << " "<< parameters.s_slave << " " << parameters.s_master <<  std::endl;

  if(dataTree == 0)
  {
    std::cout << "Tree is missing, skipping." << std::endl;
     return;
  }


  uncertain_double master_value, slave_value, residual;

  dataTree->SetBranchAddress(parameters.s_slave.c_str(), &slave );
  dataTree->SetBranchAddress(parameters.s_master.c_str(), &master );
  dataTree->SetBranchAddress("slave_cov_X_2",  &slave_cov_x);
  dataTree->SetBranchAddress("slave_cov_Y_2",  &slave_cov_y);
  dataTree->SetBranchAddress("slave_chi2",  &slave_chi2);
  dataTree->SetBranchAddress("slave_ndf",  &slave_ndf);
  if(useRecVdrift)
  {
    dataTree->SetBranchAddress("slave_recVDrift",  &slave_vdrift);
  }
  dataTree->SetBranchAddress("slave_recChamberXCenter",  &slave_chamberxcenter);
  dataTree->SetBranchAddress("slave_recChamberYCenter",  &slave_chamberycenter);
  dataTree->SetBranchAddress("slave_recChamberZCenter",  &slave_chamberzcenter);
  dataTree->SetBranchAddress("slave_recChamberRotAngleX",  &slave_chamberrotanglex);
  dataTree->SetBranchAddress("slave_recChamberRotAngleY",  &slave_chamberrotangley);
  dataTree->SetBranchAddress("slave_recChamberRotAngleZ",  &slave_chamberrotanglez);

  bool tofSource = true;
  if(parameters.s_master!="tofX" && parameters.s_master!="tofY"){
    dataTree->SetBranchAddress("master_cov_X_2",  &master_cov_x);
    dataTree->SetBranchAddress("master_cov_Y_2",  &master_cov_y);
    dataTree->SetBranchAddress("master_chi2",  &master_chi2);
    dataTree->SetBranchAddress("master_ndf",  &master_ndf);
    tofSource = false;
    if(useRecVdrift) {
      dataTree->SetBranchAddress("master_recVDrift",  &master_vdrift);
    }
    dataTree->SetBranchAddress("master_recChamberXCenter",  &master_chamberxcenter);
    dataTree->SetBranchAddress("master_recChamberYCenter",  &master_chamberycenter);
    dataTree->SetBranchAddress("master_recChamberZCenter",  &master_chamberzcenter);
    dataTree->SetBranchAddress("master_recChamberRotAngleX",  &master_chamberrotanglex);
    dataTree->SetBranchAddress("master_recChamberRotAngleY",  &master_chamberrotangley);
    dataTree->SetBranchAddress("master_recChamberRotAngleZ",  &master_chamberrotanglez);
  }


  dataTree->SetBranchAddress("eventUnixTime", &eventSecond);
  dataTree->SetBranchAddress("eventNumber", &eventNumber);
  dataTree->SetBranchAddress("runNumber", &runNumber);
  dataTree->SetBranchAddress("sectionId", &section);

  std::cout << "branchs" << std::endl;

  if(needhistos==true)
    {
      double min_resy = -20;
      double max_resy = 20;
      double min_y = -270;
      double max_y = 270;
      if(useRecVdrift) {
        max_resy = -30000;
        min_resy = +30000;
        min_y = 0, max_y = 0;
      }
      parameters.h_master = new TH1F(("h_master" + uniqueName).c_str(), "master", 500, min_y, max_y);
      parameters.h_master->GetXaxis()->SetTitle((parameters.s_master+" [cm])").c_str());
      parameters.h_slave = new TH1F(("h_slave" + uniqueName).c_str(), "slave", 100, min_y, max_y);
      parameters.h_slave->GetXaxis()->SetTitle((parameters.s_slave+" [cm])").c_str());
      parameters.h_residual = new TH1F(("h_residual" + uniqueName).c_str(), "residual",
                                       100, min_resy, max_resy);
      parameters.h_residual->GetXaxis()->SetTitle(("residual" + parameters.s_slave+ "-" + parameters.s_master + "[cm]").c_str());
      parameters.h_residual_vs_slave = new TH2F(("h_residual_vs_slave" + uniqueName).c_str(),
                                                "residual vs slave", 100, min_y, max_y, 500, min_resy, max_resy);
      parameters.h_residual_vs_slave->GetXaxis()->SetTitle((parameters.s_slave+" [cm])").c_str() );
      parameters.h_residual_vs_slave->GetYaxis()->SetTitle("residual [cm]");
    }


  long int allnentries = dataTree->GetEntries();
  long int nentries =0;

  int oldsectionnumber = 0;
  unsigned int starttime = 0, oldstarttime = 0, endtime = 0;

  std::vector<straightFitter> fitter;

  for(int i=0; i!= allnentries; i++) {
    dataTree->GetEntry(i);
    double slave_t = slave;
    double master_t = master;
    double delta_t = 0;
    if(useRecVdrift)
    {
      slave_t/=slave_vdrift;
      if(!tofSource) {
        master_t/=master_vdrift;
        delta_t = slave_t - master_t;
      }
      else {
        delta_t = master_t;
      }

    }
    else{
     delta_t = slave_t - master_t;
    }

    if(section == -1)
    {
      std::cout << "unprocessed data, sectionId = -1" << std::endl;
    }

    if(section==0 && oldsectionnumber==0 ) {
      oldsectionnumber = section;
      continue;
    }

    if(section!=0 && section!=oldsectionnumber) {
      parameters.sections.push_back(alignmentparams::alignment_section());
      parameters.sections.back().starttime = eventSecond;
      if ( parameters.runNumbers.empty() || parameters.runNumbers.back()!=runNumber ) {
        parameters.runNumbers.push_back(runNumber);
        parameters.recChamberXCenters.push_back(slave_chamberxcenter);
        parameters.recChamberYCenters.push_back(slave_chamberycenter);
        parameters.recChamberZCenters.push_back(slave_chamberzcenter);
        parameters.recChamberRotAngleXs.push_back(slave_chamberrotanglex);
        parameters.recChamberRotAngleYs.push_back(slave_chamberrotangley);
        parameters.recChamberRotAngleZs.push_back(slave_chamberrotanglez);
      }
      fitter.push_back(straightFitter());
      parameters.sections.back().n_entries=0;
      parameters.sections.back().n_dropped_entries=0;
    }

    if(section!=0) {
      /// the endtime variable updated in every entry for this section, until the last entry is filled, so at the end it will contain the ending time
      parameters.sections.back().endtime = eventSecond;
      parameters.sections.back().n_entries++;

      //should read the correlation term too
      if( TMath::Sqrt( slave_cov_x + slave_cov_y   ) > 2
     	  || slave_ndf < min_slave_ndf
        || (  !tofSource && (  master_ndf < min_master_ndf || TMath::Sqrt( master_cov_x + master_cov_y  ) > 2 ) ) )
      {
	/*        parameters.sections.back().n_dropped_entries++;
                  std::cout << " uniqueName " << uniqueName
                  << " slave_cov_x " << slave_cov_x
                  << " slave_cov_y " << slave_cov_y
                  << " slave_ndf " << slave_ndf << " " << min_slave_ndf
                  << " master_cov_x " << master_cov_x
                  << " master_cov_y " << master_cov_y
                  << " master_ndf " << master_ndf << " " << min_master_ndf << std::endl;
	*/
        nbad_tracks++;
        oldsectionnumber = section;
        continue;
      }

      fitter.back().AddEntry(slave_t, delta_t);

      parameters.sections.back().residual+=delta_t;
      parameters.sections.back().master+=master_t;
      parameters.sections.back().slave+=slave_t;

    }

    if(needhistos){
      parameters.h_slave->Fill(slave_t);
      parameters.h_master->Fill(master_t);
      parameters.h_residual->Fill(delta_t);
      parameters.h_residual_vs_slave->Fill(slave_t, delta_t);
    }

    oldsectionnumber = section;
  }

  std::cout << "dropped " << nbad_tracks << " tracks from " << allnentries << std::endl;

  unsigned int size = parameters.sections.size();
  double xmiddle[size];
  double xlow[size];
  double xhigh[size];

  double ymiddle[size];
  double ylow[size];
  double yhigh[size];

  double ymiddle2[size];
  double ylow2[size];
  double yhigh2[size];

  for(int i=0; i!= size; i++)
    {
      if(parameters.sections[i].n_entries - parameters.sections[i].n_dropped_entries > 3){
	fitter[i].GetResults(parameters.sections[i].offset,
			     parameters.sections[i].slope,
			     parameters.sections[i].cov_offset,
			     parameters.sections[i].cov_slope,
			     parameters.sections[i].cov_slope_offset);
	parameters.sections[i].Print();
      } else {
        std::cout << "Number of remaining entries in section is less then 4" << std::endl;
      }


/*      xmiddle[i] = 0.5*(parameters.sections[i].starttime+parameters.sections[i].endtime);
      xlow[i] = 0.5 * (parameters.sections[i].endtime-parameters.sections[i].starttime);
      xhigh[i] = xlow[i];

      ymiddle[i] = parameters.sections[i].slope;
      ylow[i] = sqrt(parameters.sections[i].cov_slope);
      yhigh[i] = ylow[i];

      ymiddle2[i] = parameters.sections[i].offset;
      ylow2[i] = sqrt(parameters.sections[i].cov_offset);
      yhigh2[i] = ylow2[i];*/
    }
/*
  parameters.g_slope = new TGraphAsymmErrors(size, xmiddle, ymiddle, xlow, xhigh, ylow, yhigh);
  parameters.g_slope->GetXaxis()->SetTitle("eventSecond");
  parameters.g_slope->GetYaxis()->SetTitle("fit slope");
  parameters.g_slope->SetTitle(("Fit slope for "+uniqueName
                                +" "+ parameters.s_slave+" - "
                                +parameters.s_master+" vs "
                                +parameters.s_master).c_str());
  parameters.g_offset = new TGraphAsymmErrors(size, xmiddle, ymiddle2, xlow, xhigh, ylow2, yhigh2);
  parameters.g_offset->GetXaxis()->SetTitle("eventSecond");
  parameters.g_offset->GetYaxis()->SetTitle("fit offset");
  parameters.g_offset->SetTitle(("Fit offset for "+uniqueName
                                +" "+ parameters.s_slave+" - "
                                +parameters.s_master+" vs "
                                +parameters.s_master).c_str());
*/
}

void VDriftFromTof(alignmentparams const & mainTPCparams, vdrift_corrections &mainTPCvdrift, bool useRecVdrift = true )
{
  mainTPCvdrift.sections.clear();
  for(int i=0; i!= mainTPCparams.sections.size(); i++) {
    vdrift_corrections::vdrift_section new_section;
    if((double)mainTPCparams.sections[i].n_dropped_entries/(double)mainTPCparams.sections[i].n_entries < maxInvalidEntryRate){
      if(useRecVdrift) {
        new_section.vdrift_scale = mainTPCparams.sections[i].slope;
      }
      else {
         new_section.vdrift_scale = 1.0 - mainTPCparams.sections[i].slope;
      }
      new_section.vdrift_scale_error = mainTPCparams.sections[i].cov_slope;
      new_section.starttime = mainTPCparams.sections[i].starttime;
      new_section.endtime = mainTPCparams.sections[i].endtime;


//     std::cout << "start, end = " << mainTPCparams.sections[i].starttime << " " << mainTPCparams.sections[i].endtime << std::endl;
//     std::cout << " alignment slope = " << mainTPCparams.sections[i].slope << std::endl;


      mainTPCvdrift.sections.push_back(new_section);
    }
  }
}

void Y0FromTof(alignmentparams const & mainTPCparams, y0_corrections &mainTPCy_shift )
{
  mainTPCy_shift.sections.clear();
  for(int i=0; i!= mainTPCparams.sections.size(); i++) {
    y0_corrections::y_shift_section new_section;
    if((double)mainTPCparams.sections[i].n_dropped_entries/(double)mainTPCparams.sections[i].n_entries < maxInvalidEntryRate){
      //We have to add y_shift to the used detector position, so this has opposite sign, compared to the fit offset
      new_section.y_shift = -mainTPCparams.sections[i].offset;
      new_section.y_shift_error = mainTPCparams.sections[i].cov_offset;
      new_section.starttime = mainTPCparams.sections[i].starttime;
      new_section.endtime = mainTPCparams.sections[i].endtime;
      mainTPCy_shift.sections.push_back(new_section);
    }
  }
  mainTPCy_shift.recGeomParameters.runNumberMin = (!mainTPCparams.runNumbers.empty() ? mainTPCparams.runNumbers.front() : 0);
  mainTPCy_shift.recGeomParameters.runNumberMax = (!mainTPCparams.runNumbers.empty() ? mainTPCparams.runNumbers.back() : 0);
  mainTPCy_shift.recGeomParameters.chamberXCenter = (!mainTPCparams.recChamberXCenters.empty() ? mainTPCparams.recChamberXCenters.front() : 0.0);
  mainTPCy_shift.recGeomParameters.chamberYCenter = (!mainTPCparams.recChamberYCenters.empty() ? mainTPCparams.recChamberYCenters.front() : 0.0);
  mainTPCy_shift.recGeomParameters.chamberZCenter = (!mainTPCparams.recChamberZCenters.empty() ? mainTPCparams.recChamberZCenters.front() : 0.0);
  mainTPCy_shift.recGeomParameters.chamberRotAngleX = (!mainTPCparams.recChamberRotAngleXs.empty() ? mainTPCparams.recChamberRotAngleXs.front() : 0.0);
  mainTPCy_shift.recGeomParameters.chamberRotAngleY = (!mainTPCparams.recChamberRotAngleYs.empty() ? mainTPCparams.recChamberRotAngleYs.front() : 0.0);
  mainTPCy_shift.recGeomParameters.chamberRotAngleZ = (!mainTPCparams.recChamberRotAngleZs.empty() ? mainTPCparams.recChamberRotAngleZs.front() : 0.0);
}

void VDriftFromTPC(alignmentparams const & slaveTPCmasterTPCparams,
                   vdrift_corrections const & masterTPCvdrift,
                   vdrift_corrections &slaveTPCvdrift, bool swapSlaveMaster = false)
{
  slaveTPCvdrift.sections.clear();
  for(int i=0; i!= slaveTPCmasterTPCparams.sections.size(); i++) {
    for(int j=0; j!= masterTPCvdrift.sections.size(); j++ ) {
      if(masterTPCvdrift.sections[j].endtime < slaveTPCmasterTPCparams.sections[i].starttime) { continue;}
      if(masterTPCvdrift.sections[j].starttime > slaveTPCmasterTPCparams.sections[i].endtime) { break;}
      unsigned int starttime;
      unsigned int endtime;
      if(!((double)slaveTPCmasterTPCparams.sections[i].n_dropped_entries/(double)slaveTPCmasterTPCparams.sections[i].n_entries < maxInvalidEntryRate)){
	std::cout << slaveTPCmasterTPCparams.sections[i].n_dropped_entries << " " << slaveTPCmasterTPCparams.sections[i].n_entries << std::endl;
	continue;}
      if( masterTPCvdrift.sections[j].starttime < slaveTPCmasterTPCparams.sections[i].starttime ) {
        starttime = slaveTPCmasterTPCparams.sections[i].starttime;
      }
      else {
        starttime = masterTPCvdrift.sections[j].starttime;
      }

      if(masterTPCvdrift.sections[j].endtime < slaveTPCmasterTPCparams.sections[i].endtime ) {
        endtime = masterTPCvdrift.sections[j].endtime;
      }
      else {
        endtime = slaveTPCmasterTPCparams.sections[i].endtime;
      }

      vdrift_corrections::vdrift_section new_section;
      new_section.starttime = starttime;
      new_section.endtime = endtime;

//       std::cout << "start, end = " << starttime << " " << endtime << std::endl;
//       std::cout << " alignment slope = " << slaveTPCmasterTPCparams.sections[i].slope << std::endl;
//       std::cout << " master   vdrift = " << masterTPCvdrift.sections[j].vdrift_scale << std::endl;

      if(!swapSlaveMaster) {
        double const rel_vdrift_scale = 1.0 - slaveTPCmasterTPCparams.sections[i].slope;
        new_section.vdrift_scale_error = sqr(rel_vdrift_scale)*masterTPCvdrift.sections[j].vdrift_scale_error
                                      + sqr(masterTPCvdrift.sections[j].vdrift_scale)*slaveTPCmasterTPCparams.sections[i].cov_slope;
        new_section.vdrift_scale = rel_vdrift_scale*masterTPCvdrift.sections[j].vdrift_scale;
      }
      else {
        double const rel_vdrift_scale = 1.0 - slaveTPCmasterTPCparams.sections[i].slope;
        new_section.vdrift_scale = masterTPCvdrift.sections[j].vdrift_scale/rel_vdrift_scale;
        new_section.vdrift_scale_error = sqr(masterTPCvdrift.sections[j].vdrift_scale / sqr(rel_vdrift_scale) ) * slaveTPCmasterTPCparams.sections[i].cov_slope
                                        + sqr(1.0/rel_vdrift_scale) * masterTPCvdrift.sections[j].vdrift_scale_error;
      }
      slaveTPCvdrift.sections.push_back(new_section);
    }
  }
}

void Y0FromTPC(alignmentparams const & slaveTPCmasterTPCparams,
               y0_corrections const & masterTPCy_shifts,
               y0_corrections &slaveTPCy_shifts, bool swapSlaveMaster = false)
{
  slaveTPCy_shifts.sections.clear();
  for(int i=0; i!= slaveTPCmasterTPCparams.sections.size(); i++) {
    for(int j=0; j!= masterTPCy_shifts.sections.size(); j++ ) {
      if(masterTPCy_shifts.sections[j].endtime <  slaveTPCmasterTPCparams.sections[i].starttime ) { continue;}
      if(masterTPCy_shifts.sections[j].starttime > slaveTPCmasterTPCparams.sections[i].endtime) { break;}
      unsigned int starttime;
      unsigned int endtime;
      if(!((double)slaveTPCmasterTPCparams.sections[i].n_dropped_entries/(double)slaveTPCmasterTPCparams.sections[i].n_entries < maxInvalidEntryRate)){continue;}
      if(masterTPCy_shifts.sections[j].starttime < slaveTPCmasterTPCparams.sections[i].starttime ) {
        starttime = slaveTPCmasterTPCparams.sections[i].starttime;
      }
      else {
        starttime = masterTPCy_shifts.sections[j].starttime;
      }

      if(masterTPCy_shifts.sections[j].endtime < slaveTPCmasterTPCparams.sections[i].endtime ) {
        endtime = masterTPCy_shifts.sections[j].endtime;
      }
      else {
        endtime = slaveTPCmasterTPCparams.sections[i].endtime;
      }

      y0_corrections::y_shift_section new_section;
      new_section.starttime = starttime;
      new_section.endtime = endtime;

      if(!swapSlaveMaster) {
        new_section.y_shift = masterTPCy_shifts.sections[j].y_shift - slaveTPCmasterTPCparams.sections[i].offset;
        new_section.y_shift_error = masterTPCy_shifts.sections[j].y_shift_error + slaveTPCmasterTPCparams.sections[i].cov_offset;
      }
      else {
        new_section.y_shift = -masterTPCy_shifts.sections[j].y_shift - slaveTPCmasterTPCparams.sections[i].offset;
        new_section.y_shift_error = masterTPCy_shifts.sections[j].y_shift_error + slaveTPCmasterTPCparams.sections[i].cov_offset;
      }
      slaveTPCy_shifts.sections.push_back(new_section);
    }
  }

  slaveTPCy_shifts.recGeomParameters.runNumberMin = (!slaveTPCmasterTPCparams.runNumbers.empty() ? slaveTPCmasterTPCparams.runNumbers.front() : 0);
  slaveTPCy_shifts.recGeomParameters.runNumberMax = (!slaveTPCmasterTPCparams.runNumbers.empty() ? slaveTPCmasterTPCparams.runNumbers.back() : 0);
  slaveTPCy_shifts.recGeomParameters.chamberXCenter = (!slaveTPCmasterTPCparams.recChamberXCenters.empty() ? slaveTPCmasterTPCparams.recChamberXCenters.front() : 0.0);
  slaveTPCy_shifts.recGeomParameters.chamberYCenter = (!slaveTPCmasterTPCparams.recChamberYCenters.empty() ? slaveTPCmasterTPCparams.recChamberYCenters.front() : 0.0);
  slaveTPCy_shifts.recGeomParameters.chamberZCenter = (!slaveTPCmasterTPCparams.recChamberZCenters.empty() ? slaveTPCmasterTPCparams.recChamberZCenters.front() : 0.0);
  slaveTPCy_shifts.recGeomParameters.chamberRotAngleX = (!slaveTPCmasterTPCparams.recChamberRotAngleXs.empty() ? slaveTPCmasterTPCparams.recChamberRotAngleXs.front() : 0.0);
  slaveTPCy_shifts.recGeomParameters.chamberRotAngleY = (!slaveTPCmasterTPCparams.recChamberRotAngleYs.empty() ? slaveTPCmasterTPCparams.recChamberRotAngleYs.front() : 0.0);
  slaveTPCy_shifts.recGeomParameters.chamberRotAngleZ = (!slaveTPCmasterTPCparams.recChamberRotAngleZs.empty() ? slaveTPCmasterTPCparams.recChamberRotAngleZs.front() : 0.0);
}

void PrintDriftVelocityScales(vdrift_corrections const & TPCvdrifts )
{
  std::cout << " Drift velocities scaling for " << TPCvdrifts.tpcName << std::endl;
  for(int i=0; i!= TPCvdrifts.sections.size(); i++)
  {
    std::cout << "   Time period from " << TPCvdrifts.sections[i].starttime << " to " << TPCvdrifts.sections[i].endtime << std::endl;
    std::cout << "    scaling = " << TPCvdrifts.sections[i].vdrift_scale << " +-" << sqrt( TPCvdrifts.sections[i].vdrift_scale_error ) << std::endl;
  }
}

void PrintY0Shifts(y0_corrections const & TPCy0s)
{
  std::cout << " Vertical shifts for " << TPCy0s.tpcName << std::endl;
  for(int i=0; i!= TPCy0s.sections.size(); i++)
  {
    std::cout << "   Time period from " <<  TPCy0s.sections[i].starttime << " to " << TPCy0s.sections[i].endtime << std::endl;
    std::cout << "    shift = " << TPCy0s.sections[i].y_shift << " +-" << sqrt( TPCy0s.sections[i].y_shift_error ) << std::endl;
  }
}

void PlotDriftVelocities(vdriftparams &TPCvdrifts, int timeOffset = 0)
{
  unsigned int size = TPCvdrifts.valid_size;

  double xmiddle[size];
  double xlow[size];
  double xhigh[size];

  double ymiddle[size];
  double ylow[size];
  double yhigh[size];

  unsigned int vdriftReferenceIndex[size];

  double const factor = 3600;

  std::vector<double> xv(size), yv(size), yError(size);

  int j=0;
  for(int i=0; i!= TPCvdrifts.sections.size(); i++)
  {
    if(TPCvdrifts.sections[i].validity){
      xmiddle[j] = ( 0.5*(TPCvdrifts.sections[i].starttime+TPCvdrifts.sections[i].endtime) - timeOffset )/factor;
      xlow[j] = 0.5 * (TPCvdrifts.sections[i].endtime-TPCvdrifts.sections[i].starttime)/factor;
      xhigh[j] = xlow[j];

      ymiddle[j] = TPCvdrifts.sections[i].vdrift;
      ylow[j] = sqrt(TPCvdrifts.sections[i].vdrift_error);
      yhigh[j] = ylow[j];


      xv[j] = xmiddle[j];
      yv[j] = ymiddle[j];
      yError[j] = ylow[j];

      vdriftReferenceIndex[j] = i;
      j++;
    }
  }
  TPCvdrifts.g_slope = new TGraphAsymmErrors(size, xmiddle, ymiddle, xlow, xhigh, ylow, yhigh);
  std::stringstream s_timeOffset;
  s_timeOffset << timeOffset;
  std::stringstream s_factor;
  s_factor << factor;
  if(factor==1.0)
  {
    TPCvdrifts.g_slope->GetXaxis()->SetTitle(("eventUnixTime - " + s_timeOffset.str() ).c_str() );
  }
  else
  {
    TPCvdrifts.g_slope->GetXaxis()->SetTitle(("(eventUnixTime - " + s_timeOffset.str() + ")/" + s_factor.str()).c_str() );
  }
  TPCvdrifts.g_slope->GetYaxis()->SetTitle("vdrift");

  TPCvdrifts.g_slope->SetTitle(("Drift velocity of "+TPCvdrifts.tpcName).c_str());

  if(startRange != stopRange) {
    TPCvdrifts.g_slope->GetXaxis()->SetRangeUser(startRange, stopRange);
  }


  double eventUnixTimes[TPCvdrifts.recVdrifts.size()];
  for(int i=0; i!= TPCvdrifts.recVdrifts.size(); i++)
    {
      eventUnixTimes[i] = (TPCvdrifts.eventUnixTime[i] - timeOffset)/factor;
    }
  TPCvdrifts.g_recVdrift = new TGraph(TPCvdrifts.recVdrifts.size(),
				      eventUnixTimes,
				      &TPCvdrifts.recVdrifts[0]);
  TPCvdrifts.g_recVdrift->SetLineColor(kRed);


  //fit
  fitter lineFitter(xv, yv, yError );

  //TPCvdrifts.linearFits = lineFitter.minimalExpectedChi2Fit(9);
  TPCvdrifts.linearFits = lineFitter.orderedFit(3);
  TPCvdrifts.smoothFit = new multiLineWrapper(TPCvdrifts.linearFits, xv);
  TPCvdrifts.smoothFunctor = new ROOT::Math::ParamFunctor(TPCvdrifts.smoothFit, &multiLineWrapper::eval);
  TPCvdrifts.smoothFunc = new TF1(("smootFunc"+TPCvdrifts.tpcName).c_str() , *TPCvdrifts.smoothFunctor, xv.front(), xv.back(), 0);
  for(unsigned int j=0; j!= TPCvdrifts.linearFits.size(); j++) {
    TPCvdrifts.linearFits[j].vdriftStartIndex = vdriftReferenceIndex[ TPCvdrifts.linearFits[j].points.front() ];
    TPCvdrifts.linearFits[j].vdriftStopIndex = vdriftReferenceIndex[ TPCvdrifts.linearFits[j].points.back() ];
  }
  TPCvdrifts.timeOffset = timeOffset;
  TPCvdrifts.factor = factor;
}


void PlotDriftVelocityScales(vdrift_corrections &TPCvdrifts , int timeOffset = 0 )
{
  unsigned int size = TPCvdrifts.sections.size();
  double xmiddle[size];
  double xlow[size];
  double xhigh[size];

  double ymiddle[size];
  double ylow[size];
  double yhigh[size];

  double const factor = 3600;

  for(int i=0; i!= size; i++)
  {
    xmiddle[i] = ( 0.5*(TPCvdrifts.sections[i].starttime+TPCvdrifts.sections[i].endtime) - timeOffset )/factor;
    xlow[i] = 0.5 * (TPCvdrifts.sections[i].endtime-TPCvdrifts.sections[i].starttime)/factor;
    xhigh[i] = xlow[i];

    ymiddle[i] = TPCvdrifts.sections[i].vdrift_scale;
    ylow[i] = sqrt(TPCvdrifts.sections[i].vdrift_scale_error);
    yhigh[i] = ylow[i];

  }
  TPCvdrifts.g_vdrift_correction = new TGraphAsymmErrors(size, xmiddle, ymiddle, xlow, xhigh, ylow, yhigh);
  std::stringstream s_timeOffset;
  s_timeOffset << timeOffset;
  std::stringstream s_factor;
  s_factor << factor;
  if(factor==1.0)
  {
    TPCvdrifts.g_vdrift_correction->GetXaxis()->SetTitle(("eventUnixTime - " + s_timeOffset.str() ).c_str() );
  }
  else
  {
    TPCvdrifts.g_vdrift_correction->GetXaxis()->SetTitle(("(eventUnixTime - " + s_timeOffset.str() + ")/" + s_factor.str()).c_str() );
  }
  TPCvdrifts.g_vdrift_correction->GetYaxis()->SetTitle("vdrift correction factor");

  TPCvdrifts.g_vdrift_correction->SetTitle(("Drift velocity correction factor for "+TPCvdrifts.tpcName).c_str());
  if(startRange != stopRange) {
    TPCvdrifts.g_vdrift_correction->GetXaxis()->SetRangeUser(startRange, stopRange);
  }
}

void PlotY0Shifts(y0_corrections &TPCy0s, int timeOffset)
{
  unsigned int size = TPCy0s.sections.size();
  double xmiddle[size];
  double xlow[size];
  double xhigh[size];

  double ymiddle[size];
  double ylow[size];
  double yhigh[size];

  std::vector<double> xv, yv, yError;
  std::vector<unsigned int> y0ReferenceIndex;


  double const factor = 3600;
  //  int j=0;
  std::cout << "Generating Y0 shift plot for " << TPCy0s.tpcName << std::endl;
  for(int i=0; i!= size; i++)
  {
    xmiddle[i] = ( 0.5*(TPCy0s.sections[i].starttime+TPCy0s.sections[i].endtime) - timeOffset )/factor;
    xlow[i] = 0.5 * (TPCy0s.sections[i].endtime-TPCy0s.sections[i].starttime)/factor;
    xhigh[i] = xlow[i];

    ymiddle[i] = TPCy0s.sections[i].y_shift;
    ylow[i] = sqrt(TPCy0s.sections[i].y_shift_error);
    yhigh[i] = ylow[i];

    if(TMath::Finite(ymiddle[i]) && TMath::Finite(ylow[i]) && TMath::Finite(xmiddle[i])) {
      xv.push_back(xmiddle[i]);
      yv.push_back(ymiddle[i]);
      yError.push_back(ylow[i]);
      y0ReferenceIndex.push_back(i); //in case a section is invalid, j!=i . Need to check if section was valid...
      //j++;
    }
  }
  //std::cout << "finished illing arrays for " << TPCy0s.tpcName << std::endl;
  TPCy0s.g_y0_correction= new TGraphAsymmErrors(size, xmiddle, ymiddle, xlow, xhigh, ylow, yhigh);
  std::stringstream s_timeOffset;
  s_timeOffset << timeOffset;
  std::stringstream s_factor;
  s_factor << factor;
  if(factor==1.0)
  {
    TPCy0s.g_y0_correction->GetXaxis()->SetTitle(("eventUnixTime - " + s_timeOffset.str() ).c_str() );
  }
  else
  {
    TPCy0s.g_y0_correction->GetXaxis()->SetTitle(("(eventUnixTime - " + s_timeOffset.str() + ")/" + s_factor.str()).c_str() );
  }
  TPCy0s.g_y0_correction->GetYaxis()->SetTitle("vertical shift [cm]");

  TPCy0s.g_y0_correction->SetTitle(("Vertical shift for "+TPCy0s.tpcName).c_str());
  if(startRange != stopRange) {
    TPCy0s.g_y0_correction->GetXaxis()->SetRangeUser(startRange, stopRange);
  }

  if(xv.size()>1) {
    fitter lineFitter(xv, yv, yError);
    TPCy0s.linearFits = lineFitter.orderedFit(3);
    TPCy0s.smoothFit = new multiLineWrapper(TPCy0s.linearFits, xv);
    TPCy0s.smoothFunctor = new ROOT::Math::ParamFunctor(TPCy0s.smoothFit, &multiLineWrapper::eval);
    TPCy0s.smoothFunc = new TF1(("smootFunc"+TPCy0s.tpcName).c_str() , *TPCy0s.smoothFunctor, xv.front(), xv.back(), 0);
    for(unsigned int j=0; j!= TPCy0s.linearFits.size(); j++) {
      TPCy0s.linearFits[j].vdriftStartIndex = y0ReferenceIndex[ TPCy0s.linearFits[j].points.front() ];
      TPCy0s.linearFits[j].vdriftStopIndex = y0ReferenceIndex[ TPCy0s.linearFits[j].points.back() ];
    }
    TPCy0s.smoothFunc->SetLineColor(kGreen);
    TPCy0s.smoothFunc->SetLineStyle(2);
  }
  TPCy0s.timeOffset = timeOffset;
  TPCy0s.factor = factor;
}

void PlotRelativeAlignment(alignmentparams &relativeAlignment, unsigned int timeOffset )
{
  unsigned int size = relativeAlignment.sections.size();
  double xmiddle[size];
  double xlow[size];
  double xhigh[size];

  double ymiddle[size];
  double ylow[size];
  double yhigh[size];

  double ymiddle2[size];
  double ylow2[size];
  double yhigh2[size];
  double const factor = 3600;
  for(int i=0; i!= size; i++)
  {
    xmiddle[i] = (0.5*(relativeAlignment.sections[i].starttime+relativeAlignment.sections[i].endtime) - timeOffset)/factor ;
    xlow[i] = 0.5 * (relativeAlignment.sections[i].endtime-relativeAlignment.sections[i].starttime)/factor;
    xhigh[i] = xlow[i];

    ymiddle[i] = relativeAlignment.sections[i].slope;
    ylow[i] = sqrt(relativeAlignment.sections[i].cov_slope);
    yhigh[i] = ylow[i];

    ymiddle2[i] = relativeAlignment.sections[i].offset;
    ylow2[i] = sqrt(relativeAlignment.sections[i].cov_offset);
    yhigh2[i] = ylow2[i];
  }
  std::stringstream s_timeOffset;
  s_timeOffset << timeOffset;
  std::stringstream s_factor;
  s_factor << factor;

  relativeAlignment.g_slope = new TGraphAsymmErrors(size, xmiddle, ymiddle, xlow, xhigh, ylow, yhigh);
  if(factor==1.0)
  {
    relativeAlignment.g_slope->GetXaxis()->SetTitle(("eventUnixTime - " + s_timeOffset.str() ).c_str() );
  }
  else
  {
    relativeAlignment.g_slope->GetXaxis()->SetTitle(("(eventUnixTime - " + s_timeOffset.str() + ")/" + s_factor.str()).c_str() );
  }

  relativeAlignment.g_slope->GetYaxis()->SetTitle("fit slope");
  relativeAlignment.g_slope->SetTitle(("Fit slope for "+ relativeAlignment.s_TPCvsTPC
                                +" "+ relativeAlignment.s_slave+" - "
                                +relativeAlignment.s_master+" vs "
                                +relativeAlignment.s_master).c_str());
  relativeAlignment.g_offset = new TGraphAsymmErrors(size, xmiddle, ymiddle2, xlow, xhigh, ylow2, yhigh2);

  if(factor==1.0)
  {
    relativeAlignment.g_offset->GetXaxis()->SetTitle(("eventUnixTime - " + s_timeOffset.str() ).c_str() );
  }
  else
  {
    relativeAlignment.g_offset->GetXaxis()->SetTitle(("(eventUnixTime - " + s_timeOffset.str() + ")/" + s_factor.str()).c_str() );
  }

  relativeAlignment.g_offset->GetYaxis()->SetTitle("fit offset");
  relativeAlignment.g_offset->SetTitle(("Fit offset for "+relativeAlignment.s_TPCvsTPC
                                +" "+ relativeAlignment.s_slave+" - "
                                +relativeAlignment.s_master+" vs "
                                +relativeAlignment.s_master).c_str());

  if(startRange!=stopRange) {
    relativeAlignment.g_slope->GetXaxis()->SetRangeUser(startRange, stopRange);
    relativeAlignment.g_offset->GetXaxis()->SetRangeUser(startRange, stopRange);
  }
}

void PrettyPlotAlignmentHisto(alignmentparams const & relativeAlignment)
{
  std::string canvasName_ = "c_" + relativeAlignment.s_TPCvsTPC + "_histos";
  std::string canvasName = relativeAlignment.s_TPCvsTPC + " histos";
  TCanvas * temp = new TCanvas(canvasName_.c_str(), canvasName.c_str(),1200, 800 );
  temp->Divide(2,2);
  temp->cd(1);
  if(relativeAlignment.h_master!=0) relativeAlignment.h_master->Draw();
  temp->cd(2);
  if(relativeAlignment.h_slave!=0) relativeAlignment.h_slave->Draw();
  temp->cd(3);
  if(relativeAlignment.h_residual!=0) relativeAlignment.h_residual->Draw();
  temp->cd(4);
  if(relativeAlignment.h_residual_vs_slave!=0) relativeAlignment.h_residual_vs_slave->Draw();



}

void PrettyPlot(TGraphAsymmErrors  *mtpcl,
                TGraphAsymmErrors  *mtpcr,
                TGraphAsymmErrors  *vtpc1,
                TGraphAsymmErrors  *vtpc2,
                TGraphAsymmErrors  *gtpc,
                TGraphAsymmErrors  *mtpcr_crosscheck,
                TGraphAsymmErrors  *gtpc_crosscheck,
                std::string const &calibValueName,
                std::string const &calibValuePrintName,
		TF1 * mtpcl_func = 0,
		TF1 * mtpcr_func = 0,
		TF1 * vtpc1_func = 0,
		TF1 * vtpc2_func = 0,
		TF1 * gtpc_func = 0,
		TF1 * mtpcr_crosscheck_func = 0,
		TF1 * gtpc_crosscheck_func = 0
               )
{
  std::string canvasName = "c_"+calibValuePrintName+ "_ML_V1_V2_G_pretty";
  TCanvas * temp1 = new TCanvas(canvasName.c_str(),
                                ("MTPCL VTPC1 VTPC2 GTPC "+calibValueName).c_str(), 1200, 800);
  temp1->Divide(2,2);
  temp1->cd(1);
  if(mtpcl!=0) {
    mtpcl->Draw("ap");
    if(mtpcl_func!=0) mtpcl_func->Draw("l same");
  }
  temp1->cd(2);
  if(vtpc2!=0) {
    vtpc2->Draw("ap");
    if(vtpc2_func!=0) vtpc2_func->Draw("l same");
  }
  temp1->cd(3);
  if(vtpc1!=0) {
    vtpc1->Draw("ap");
    if(vtpc1_func!=0) vtpc1_func->Draw("l same");
  }
  temp1->cd(4);
  if(gtpc!=0) {
    gtpc->Draw("ap");
    if(gtpc_func!=0) gtpc_func->Draw("l same");
  }

  temp1->Print((canvasName + ".eps").c_str());

  std::string canvasName2 = "c_" + calibValuePrintName + "_MR_MRcrosscheck_Gcrosscheck_pretty";


  TCanvas * temp2 = new TCanvas(canvasName2.c_str(),
                                ("MTPCR MTPCR_crosscheck GTPC_crosscheck"+calibValueName).c_str(), 1200, 800);

  temp2->Divide(2,2);

  temp2->cd(1);
  if(mtpcr!=0) {
    mtpcr->Draw("ap");
    if(mtpcr_func!=0) mtpcr_func->Draw("l same");
  }
  temp2->cd(2);
  if(mtpcr_crosscheck!=0) {
    mtpcr_crosscheck->Draw("ap");
    if(mtpcr_crosscheck_func!=0) mtpcr_crosscheck_func->Draw("l same");
  }
  temp2->cd(3);
  if(gtpc_crosscheck!=0) {
    gtpc_crosscheck->Draw("ap");
    if(gtpc_crosscheck_func!=0) gtpc_crosscheck_func->Draw("l same");
  }


  temp2->Print((canvasName2 + ".eps").c_str());

}

void WriteAbsDriftVelocities(vdriftparams const &TPCVdrifts)
{
  std::cout << "Writing file " << TPCVdrifts.tpcName << ".txt" << std::endl;
  std::ofstream ofile( (TPCVdrifts.tpcName+".txt").c_str(), std::ofstream::out);
  ofile << "# unixTime[sec] vDrift[cm/usec]\n";
  bool writeStartValue = true;
  for(unsigned int i=0; i < TPCVdrifts.linearFits.size(); i++) {
    unsigned int const startTime = TPCVdrifts.sections[ TPCVdrifts.linearFits[i].vdriftStartIndex ].starttime;
    unsigned int const stopTime = TPCVdrifts.sections[ TPCVdrifts.linearFits[i].vdriftStopIndex ].endtime;
    double const startVdrift = TPCVdrifts.linearFits[i].f[0] + TPCVdrifts.linearFits[i].f[1]*((long int)startTime-(long int)TPCVdrifts.timeOffset)/TPCVdrifts.factor;
    double const stopVdrift = TPCVdrifts.linearFits[i].f[0] + TPCVdrifts.linearFits[i].f[1]*((long int)stopTime-(long int)TPCVdrifts.timeOffset)/TPCVdrifts.factor;
    if ( !writeStartValue )
      writeStartValue = true;
    else
      ofile << startTime << " " << startVdrift*1000.0 << std::endl;
    if ( i+1<TPCVdrifts.linearFits.size() ) {
      unsigned int const startTimeNext = TPCVdrifts.sections[ TPCVdrifts.linearFits[i+1].vdriftStartIndex ].starttime;
      unsigned int const stopTimeNext = TPCVdrifts.sections[ TPCVdrifts.linearFits[i+1].vdriftStopIndex ].endtime;
      double const startVdriftNext = TPCVdrifts.linearFits[i+1].f[0] + TPCVdrifts.linearFits[i+1].f[1]*((long int)startTime-(long int)TPCVdrifts.timeOffset)/TPCVdrifts.factor;
      double const stopVdriftNext = TPCVdrifts.linearFits[i+1].f[0] + TPCVdrifts.linearFits[i+1].f[1]*((long int)stopTime-(long int)TPCVdrifts.timeOffset)/TPCVdrifts.factor;
      if ( stopTime == startTimeNext ) {
        ofile << stopTime << " " << 0.5*(stopVdrift+startVdriftNext)*1000.0 << std::endl;
        writeStartValue = false;
      }
    }
    if ( writeStartValue )
      ofile << stopTime << " " << stopVdrift*1000.0 << std::endl;
  }
  ofile.close();
}

void WriteY0Alignment(y0_corrections const &TPCy0s)
{
  std::cout << "Writing file " << TPCy0s.tpcName << "_y0.txt" << std::endl;
  std::ofstream ofile( (TPCy0s.tpcName+"_y0.txt").c_str(), std::ofstream::out);
  ofile << "# unixTime[sec] y0Correction[cm] (y0=" << TPCy0s.recGeomParameters.chamberYCenter << "cm)\n";
  bool writeStartValue = true;
  for(unsigned int i=0; i < TPCy0s.linearFits.size(); i++) {
    unsigned int const startTime = TPCy0s.sections[ TPCy0s.linearFits[i].vdriftStartIndex ].starttime;
    unsigned int const stopTime = TPCy0s.sections[ TPCy0s.linearFits[i].vdriftStopIndex ].endtime;
    double const startY0 = TPCy0s.linearFits[i].f[0] + TPCy0s.linearFits[i].f[1]*((long int)startTime-(long int)TPCy0s.timeOffset)/TPCy0s.factor;
    double const stopY0 = TPCy0s.linearFits[i].f[0] + TPCy0s.linearFits[i].f[1]*((long int)stopTime-(long int)TPCy0s.timeOffset)/TPCy0s.factor;
    if ( !writeStartValue )
      writeStartValue = true;
    else
      ofile << startTime << " " << startY0 << std::endl;
    if ( i+1 < TPCy0s.linearFits.size() ) {
      unsigned int const startTimeNext = TPCy0s.sections[ TPCy0s.linearFits[i+1].vdriftStartIndex ].starttime;
      unsigned int const stopTimeNext = TPCy0s.sections[ TPCy0s.linearFits[i+1].vdriftStopIndex ].endtime;
      double const startY0Next = TPCy0s.linearFits[i+1].f[0] + TPCy0s.linearFits[i+1].f[1]*((long int)startTime-(long int)TPCy0s.timeOffset)/TPCy0s.factor;
      double const stopY0Next = TPCy0s.linearFits[i+1].f[0] + TPCy0s.linearFits[i+1].f[1]*((long int)stopTime-(long int)TPCy0s.timeOffset)/TPCy0s.factor;
      if ( stopTime == startTimeNext ) {
        ofile << stopTime << " " << 0.5*(stopY0+startY0Next) << std::endl;
        writeStartValue = false;
      }
    }
    if ( writeStartValue )
      ofile << stopTime << " " << stopY0 << std::endl;
  }
  ofile.close();
}


void WriteAbsY0AlignmentXML(const std::string& fileName, const std::map<eCommonPlanes, y0_corrections>& allTPCy0s)
{
  std::cout << "Writing file " << fileName << std::endl;
  std::ofstream ofile(fileName.c_str(), std::ofstream::out);
  ofile << "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n\n";
  ofile << "<TPCAlignmentGeometry\n"
        << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
        << "  xsi:noNamespaceSchemaLocation=\"[SCHEMAPATH]/TPCAlignmentGeometry_DBFormat.xsd\"\n"
        << "  positionUnit=\"cm\"\n"
        << "  angleUnit=\"deg\"\n"
        << "  rotAngleOrder=\"y-z-x\"\n"
        << "  coordinateSystem=\"global\">\n\n";
  unsigned int runNumberMin = 0;
  unsigned int runNumberMax = 0;
  for(std::map<eCommonPlanes, y0_corrections >::const_iterator it = allTPCy0s.begin(); it!= allTPCy0s.end(); ++it) {
    const y0_corrections& y0corr = it->second;
    const std::string& tpcName = y0corr.tpcName;
    if ( tpcName.find("_crosscheck") != std::string::npos )
      continue;
    if ( !runNumberMin )
      runNumberMin = y0corr.recGeomParameters.runNumberMin;
    if ( !runNumberMax )
      runNumberMax = y0corr.recGeomParameters.runNumberMax;
  }
  ofile << "  <GeometryEntry runStart=\"" << runNumberMin << "\" runStop=\"" << runNumberMax << "\">\n\n";
  for(std::map<eCommonPlanes, y0_corrections >::const_iterator it = allTPCy0s.begin(); it!= allTPCy0s.end(); ++it) {
    const y0_corrections& y0corr = it->second;
    const std::string& tpcName = y0corr.tpcName;
    if ( tpcName.find("_crosscheck") != std::string::npos )
      continue;
    ofile << "    <TPC name=\"" << tpcName << "\">\n";
    double y0Shift = 0.0;
    int nY0ShiftEntries = 0;
    for ( unsigned int i=0; i < y0corr.linearFits.size(); ++i ) {
      unsigned int const startTime = y0corr.sections[ y0corr.linearFits[i].vdriftStartIndex ].starttime;
      unsigned int const stopTime = y0corr.sections[ y0corr.linearFits[i].vdriftStopIndex ].endtime;
      double const startY0 = y0corr.linearFits[i].f[0] + y0corr.linearFits[i].f[1]*((long int)startTime-(long int)y0corr.timeOffset)/y0corr.factor;
      double const stopY0 = y0corr.linearFits[i].f[0] + y0corr.linearFits[i].f[1]*((long int)stopTime-(long int)y0corr.timeOffset)/y0corr.factor;
      y0Shift += startY0;
      ++nY0ShiftEntries;
      y0Shift += stopY0;
      ++nY0ShiftEntries;
    }
    y0Shift = (nY0ShiftEntries ? y0Shift/nY0ShiftEntries : 0.0);
    const double x0 = y0corr.recGeomParameters.chamberXCenter;
    const double y0 = y0corr.recGeomParameters.chamberYCenter;
    const double z0 = y0corr.recGeomParameters.chamberZCenter;
    const double rotanglex = y0corr.recGeomParameters.chamberRotAngleX;
    const double rotangley = y0corr.recGeomParameters.chamberRotAngleY;
    const double rotanglez = y0corr.recGeomParameters.chamberRotAngleZ;
    ofile << "      <centerX> " << x0 << " </centerX>\n";
    ofile << "      <centerY> " << y0+y0Shift << " </centerY>\n";
    ofile << "      <centerZ> " << z0 << " </centerZ>\n";
    ofile << "      <rotAngleX> " << rotanglex << " </rotAngleX>\n";
    ofile << "      <rotAngleY> " << rotangley << " </rotAngleY>\n";
    ofile << "      <rotAngleZ> " << rotanglez << " </rotAngleZ>\n";
    ofile << "    </TPC>\n\n";
  }
  ofile << "  </GeometryEntry>\n\n";
  ofile << "</TPCAlignmentGeometry>\n\n";
  ofile.close();
}


void AbsDriftVelocityCampaign(std::string const & filename)
{
  TFile * file = TFile::Open(filename.c_str(), "ro");
  std::map<std::string, vdriftparams> vdrifts;
  vdrifts["MTPCL"].tpcName = "MTPCL";
  vdrifts["MTPCR"].tpcName = "MTPCR";
  vdrifts["VTPC1"].tpcName = "VTPC1" ;
  vdrifts["VTPC2"].tpcName = "VTPC2" ;
  vdrifts["GTPC"].tpcName = "GTPC" ;
  vdriftparams tofl;
  vdriftparams tofr;
  tofl.tpcName = "TOFL";
  tofr.tpcName = "TOFR";

  std::map<eCommonPlanes, std::string > treeNames;
  treeNames[eMTPCLvsTOFL] = "MTPCLvsTOFL";
  treeNames[eMTPCRvsTOFR] = "MTPCRvsTOFR";
  treeNames[eVTPC2vsMTPCL] = "VTPC2vsMTPCL";
  treeNames[eVTPC2vsMTPCR] = "VTPC2vsMTPCR";
  treeNames[eVTPC1vsVTPC2] = "VTPC1vsVTPC2";
  treeNames[eVTPC1vsGTPC] = "VTPC1vsGTPC";
  treeNames[eGTPCvsVTPC2] = "GTPCvsVTPC2";

  std::map<eCommonPlanes, TTree *> dataTrees;
  dataTrees[eMTPCLvsTOFL]  = (TTree *)file->Get(treeNames[eMTPCLvsTOFL].c_str());
  dataTrees[eMTPCRvsTOFR]  = (TTree *)file->Get(treeNames[eMTPCRvsTOFR].c_str());
  dataTrees[eVTPC1vsVTPC2] = (TTree *)file->Get(treeNames[eVTPC1vsVTPC2].c_str());
  dataTrees[eVTPC2vsMTPCL] = (TTree *)file->Get(treeNames[eVTPC2vsMTPCL].c_str());
  dataTrees[eVTPC2vsMTPCR] = (TTree *)file->Get(treeNames[eVTPC2vsMTPCR].c_str());
  dataTrees[eVTPC1vsGTPC]  = (TTree *)file->Get(treeNames[eVTPC1vsGTPC].c_str());
  dataTrees[eGTPCvsVTPC2]  = (TTree *)file->Get(treeNames[eGTPCvsVTPC2].c_str());

  fitAbsVdrift(dataTrees[eMTPCLvsTOFL],  tofl, vdrifts["MTPCL"], false);
  fitAbsVdrift(dataTrees[eMTPCRvsTOFR],  tofr, vdrifts["MTPCR"], false);
  fitAbsVdrift(dataTrees[eVTPC2vsMTPCL], vdrifts["MTPCL"],vdrifts["VTPC2"], false);
  fitAbsVdrift(dataTrees[eVTPC1vsVTPC2], vdrifts["VTPC2"],vdrifts["VTPC1"], false);
  fitAbsVdrift(dataTrees[eVTPC1vsGTPC], vdrifts["VTPC1"],vdrifts["GTPC"], true, 40, 5);


  int timeOffset = 0;
  if(vdrifts["VTPC1"].sections.size()!=0 )
    {
      timeOffset = vdrifts["MTPCL"].sections[0].starttime;
    }


  for(std::map<std::string, vdriftparams>::iterator it = vdrifts.begin();
      it!= vdrifts.end(); it++)
    {
      if(it->second.valid_size != 0)
	{
	  PlotDriftVelocities(it->second, timeOffset);
	  TCanvas * temp = new TCanvas(("c_vdrift_"+ it->second.tpcName).c_str(),(it->second.tpcName + " drift velocity").c_str(), 1 );
	  std::cout << "TCanvas " << temp->GetName() << " with title: " << temp->GetTitle() << " was initialized." << std::endl;
	  it->second.g_slope->Draw("ap");
	  it->second.g_recVdrift->Draw("lp");
	  it->second.smoothFunc->SetLineColor(kGreen);
	  it->second.smoothFunc->SetLineStyle(2);
	  it->second.smoothFunc->Draw("l same");

	  //temp->Update(); //Leads to segmentation break
	  //temp->Print(("c_abs_vdrift_"+ it->second.tpcName + ".eps").c_str()); //Leads to segmentation break
	  //temp->SaveAs(("c_abs_vdrift_"+ it->second.tpcName + ".C").c_str());//Leads to segmentation break
	  //      TCanvas * temp2 =
	  //new TCanvas(("c_recvdrift_"+ it->second.tpcName).c_str(),(it->second.tpcName + " drift velocity").c_str(), 1 );
	}
      else
	{
	  std::cout << "No entries for " << it->second.tpcName << ". No plot will be created."<< std::endl;
	}

    }

  for(std::map<std::string, vdriftparams>::iterator it = vdrifts.begin();
      it!= vdrifts.end(); it++)
    {
      if(it->second.valid_size != 0)
	{
	  WriteAbsDriftVelocities(it->second);
	}
      else
	{
	  std::cout << "No entries for " << it->second.tpcName << ". No file will be created."<< std::endl;
	}
    }


}

void DriftVelocityY0CalibrationCampaign(std::string const & filename,
                                        bool calibrateDriftVelocity = true,
                                        bool calibrateY0 = true,
					bool useRecVdrift = true,
                                        bool prettyPlot =true)
{
  //This calibration function assumes perfect alignment (except vertical)
  TFile * file = TFile::Open(filename.c_str(), "ro");




  std::string const tofHitParameter = "tofY";
  std::string const tpcTrackPamater_slave = "slave_Y";
  std::string const tpcTrackPamater_master = "master_Y";


  std::map<eCommonPlanes, std::string > treeNames;
  treeNames[eMTPCLvsTOFL] = "MTPCLvsTOFL";
  treeNames[eMTPCRvsTOFR] = "MTPCRvsTOFR";
  treeNames[eVTPC2vsMTPCL] = "VTPC2vsMTPCL";
  treeNames[eVTPC2vsMTPCR] = "VTPC2vsMTPCR";
  treeNames[eVTPC1vsVTPC2] = "VTPC1vsVTPC2";
  treeNames[eVTPC1vsGTPC] = "VTPC1vsGTPC";
  treeNames[eGTPCvsVTPC2] = "GTPCvsVTPC2";

  std::map<eCommonPlanes, alignmentparams> verticalAlignment;
  verticalAlignment[eMTPCLvsTOFL].s_master = tofHitParameter;
  verticalAlignment[eMTPCLvsTOFL].s_slave  = tpcTrackPamater_slave;
  verticalAlignment[eMTPCRvsTOFR].s_master = tofHitParameter;
  verticalAlignment[eMTPCRvsTOFR].s_slave  = tpcTrackPamater_slave;
  verticalAlignment[eVTPC1vsVTPC2].s_master = tpcTrackPamater_master;
  verticalAlignment[eVTPC1vsVTPC2].s_slave  = tpcTrackPamater_slave;
  verticalAlignment[eVTPC2vsMTPCL].s_master = tpcTrackPamater_master;
  verticalAlignment[eVTPC2vsMTPCL].s_slave  = tpcTrackPamater_slave;
  verticalAlignment[eVTPC2vsMTPCR].s_master = tpcTrackPamater_master;
  verticalAlignment[eVTPC2vsMTPCR].s_slave  = tpcTrackPamater_slave;
  verticalAlignment[eVTPC1vsGTPC].s_master = tpcTrackPamater_master;
  verticalAlignment[eVTPC1vsGTPC].s_slave  = tpcTrackPamater_slave;
  verticalAlignment[eGTPCvsVTPC2].s_master = tpcTrackPamater_master;
  verticalAlignment[eGTPCvsVTPC2].s_slave  = tpcTrackPamater_slave;

  std::map<eCommonPlanes, TTree *> dataTrees;
  dataTrees[eMTPCLvsTOFL]  = (TTree *)file->Get(treeNames[eMTPCLvsTOFL].c_str());
  dataTrees[eMTPCRvsTOFR]  = (TTree *)file->Get(treeNames[eMTPCRvsTOFR].c_str());
  dataTrees[eVTPC1vsVTPC2] = (TTree *)file->Get(treeNames[eVTPC1vsVTPC2].c_str());
  dataTrees[eVTPC2vsMTPCL] = (TTree *)file->Get(treeNames[eVTPC2vsMTPCL].c_str());
  dataTrees[eVTPC2vsMTPCR] = (TTree *)file->Get(treeNames[eVTPC2vsMTPCR].c_str());
  dataTrees[eVTPC1vsGTPC]  = (TTree *)file->Get(treeNames[eVTPC1vsGTPC].c_str());
  dataTrees[eGTPCvsVTPC2]  = (TTree *)file->Get(treeNames[eGTPCvsVTPC2].c_str());

  bool makeHistograms = true;

  fitVdrift(dataTrees[eMTPCLvsTOFL], verticalAlignment[eMTPCLvsTOFL], makeHistograms, treeNames[eMTPCLvsTOFL], useRecVdrift );
  fitVdrift(dataTrees[eMTPCRvsTOFR], verticalAlignment[eMTPCRvsTOFR], makeHistograms, treeNames[eMTPCRvsTOFR] , useRecVdrift);
  fitVdrift(dataTrees[eVTPC1vsVTPC2], verticalAlignment[eVTPC1vsVTPC2], makeHistograms, treeNames[eVTPC1vsVTPC2], useRecVdrift);
  fitVdrift(dataTrees[eVTPC2vsMTPCL], verticalAlignment[eVTPC2vsMTPCL], makeHistograms, treeNames[eVTPC2vsMTPCL], useRecVdrift);
  fitVdrift(dataTrees[eVTPC2vsMTPCR], verticalAlignment[eVTPC2vsMTPCR], makeHistograms, treeNames[eVTPC2vsMTPCR], useRecVdrift);
  fitVdrift(dataTrees[eVTPC1vsGTPC], verticalAlignment[eVTPC1vsGTPC], makeHistograms, treeNames[eVTPC1vsGTPC], useRecVdrift , 5, 40);
  fitVdrift(dataTrees[eGTPCvsVTPC2], verticalAlignment[eGTPCvsVTPC2], makeHistograms, treeNames[eGTPCvsVTPC2], useRecVdrift , 40, 5);

  PrettyPlotAlignmentHisto( verticalAlignment[eMTPCLvsTOFL] ) ;
  PrettyPlotAlignmentHisto( verticalAlignment[eMTPCRvsTOFR] ) ;
  PrettyPlotAlignmentHisto( verticalAlignment[eVTPC2vsMTPCL] ) ;
  PrettyPlotAlignmentHisto( verticalAlignment[eVTPC2vsMTPCR] ) ;
  PrettyPlotAlignmentHisto( verticalAlignment[eVTPC1vsVTPC2] ) ;
  PrettyPlotAlignmentHisto( verticalAlignment[eVTPC1vsGTPC] ) ;
  PrettyPlotAlignmentHisto( verticalAlignment[eGTPCvsVTPC2] ) ;



  int timeOffset = 0;
    if(verticalAlignment[eVTPC1vsVTPC2].sections.size()!=0 )
    {
      timeOffset = verticalAlignment[eVTPC1vsVTPC2].sections[0].starttime;
    }


  if(calibrateDriftVelocity) {
    std::map<eCommonPlanes, vdrift_corrections > vdriftScales;
    //MTPCL driftvelocity
    VDriftFromTof( verticalAlignment[eMTPCLvsTOFL], vdriftScales[eMTPCLvsTOFL], useRecVdrift  );
    vdriftScales[eMTPCLvsTOFL].tpcName = "MTPCL";
    //MTPCR driftvelocity
    VDriftFromTof( verticalAlignment[eMTPCRvsTOFR], vdriftScales[eMTPCRvsTOFR], useRecVdrift  );
    vdriftScales[eMTPCRvsTOFR].tpcName = "MTPCR";

    //VTCP2 driftvelocity
    VDriftFromTPC( verticalAlignment[eVTPC2vsMTPCL], vdriftScales[eMTPCLvsTOFL],  vdriftScales[eVTPC2vsMTPCL] );
    vdriftScales[eVTPC2vsMTPCL].tpcName = "VTPC2";
    //VTPC1 driftvelocity
    VDriftFromTPC( verticalAlignment[eVTPC1vsVTPC2], vdriftScales[eVTPC2vsMTPCL], vdriftScales[eVTPC1vsVTPC2] );
    vdriftScales[eVTPC1vsVTPC2].tpcName = "VTPC1";
    //GTPC driftvelocity
    VDriftFromTPC( verticalAlignment[eVTPC1vsGTPC],  vdriftScales[eVTPC1vsVTPC2], vdriftScales[eVTPC1vsGTPC], true );
    vdriftScales[eVTPC1vsGTPC].tpcName = "GTPC";
    //this is to crosscheck the MTPCR drift velocities
    VDriftFromTPC( verticalAlignment[eVTPC2vsMTPCR], vdriftScales[eVTPC2vsMTPCL], vdriftScales[eMTPCRvsVTPC2_crosscheck], true );
    vdriftScales[eMTPCRvsVTPC2_crosscheck].tpcName = "MTPCR_crosscheck";
    //this is to crosscheck the GTPC drift velocities
    VDriftFromTPC( verticalAlignment[eVTPC2vsMTPCR], vdriftScales[eVTPC2vsMTPCL], vdriftScales[eGTPCvsVTPC2]);
    vdriftScales[eGTPCvsVTPC2].tpcName = "GTPC_crosscheck";



    PrintDriftVelocityScales( vdriftScales[eMTPCLvsTOFL] );
    PrintDriftVelocityScales( vdriftScales[eMTPCRvsTOFR] );
    PrintDriftVelocityScales( vdriftScales[eVTPC2vsMTPCL] );
    PrintDriftVelocityScales( vdriftScales[eVTPC1vsVTPC2] );
    PrintDriftVelocityScales( vdriftScales[eVTPC1vsGTPC] );
    PrintDriftVelocityScales( vdriftScales[eMTPCRvsVTPC2_crosscheck] );
    PrintDriftVelocityScales( vdriftScales[eGTPCvsVTPC2] );

//     int timeOffset = 0;
    if(vdriftScales[eVTPC1vsVTPC2].sections.size()!=0 )
     {
       timeOffset = vdriftScales[eVTPC1vsVTPC2].sections[0].starttime;
     }


    for(std::map<eCommonPlanes, vdrift_corrections >::iterator it = vdriftScales.begin();
        it!= vdriftScales.end(); it++)
    {

      PlotDriftVelocityScales( it->second , timeOffset);
      if(!prettyPlot) {
        TCanvas * temp = new TCanvas(("c_vdrift_scale_"+ it->second.tpcName).c_str(),(it->second.tpcName + " drift velocity").c_str(), 1 );
        it->second.g_vdrift_correction->Draw("ap");
        temp->Print(("c_vdrift_scale_"+ it->second.tpcName + ".eps").c_str());
      }
    }

    if(prettyPlot) {
      PrettyPlot(vdriftScales[eMTPCLvsTOFL].g_vdrift_correction,
               vdriftScales[eMTPCRvsTOFR].g_vdrift_correction,
               vdriftScales[eVTPC2vsMTPCL].g_vdrift_correction,
               vdriftScales[eVTPC1vsVTPC2].g_vdrift_correction,
               vdriftScales[eVTPC1vsGTPC].g_vdrift_correction,
               vdriftScales[eMTPCRvsVTPC2_crosscheck].g_vdrift_correction,
               vdriftScales[eGTPCvsVTPC2].g_vdrift_correction,
               "drift velocity factor",
               "vdrift_scale");
    }
  }

  if(calibrateY0) {
    std::map<eCommonPlanes,y0_corrections > y0shifts;
    //MTPCL y shift
    Y0FromTof( verticalAlignment[eMTPCLvsTOFL], y0shifts[eMTPCLvsTOFL] );
    y0shifts[eMTPCLvsTOFL].tpcName = "MTPCL";
    //MTPCR y shift
    Y0FromTof( verticalAlignment[eMTPCRvsTOFR], y0shifts[eMTPCRvsTOFR] );
    y0shifts[eMTPCRvsTOFR].tpcName = "MTPCR";

    //VTPC2 y shift
    Y0FromTPC( verticalAlignment[eVTPC2vsMTPCL], y0shifts[eMTPCLvsTOFL], y0shifts[eVTPC2vsMTPCL] );
    y0shifts[eVTPC2vsMTPCL].tpcName = "VTPC2";
    //VTPC1 y shift
    Y0FromTPC( verticalAlignment[eVTPC1vsVTPC2], y0shifts[eVTPC2vsMTPCL], y0shifts[eVTPC1vsVTPC2] );
    y0shifts[eVTPC1vsVTPC2].tpcName = "VTPC1";
    //GTPC y shift
    Y0FromTPC( verticalAlignment[eVTPC1vsGTPC], y0shifts[eVTPC1vsVTPC2], y0shifts[eVTPC1vsGTPC] , true );
    y0shifts[eVTPC1vsGTPC].tpcName = "GTPC";
    //MTPCR y shift crosscheck from VTPC2
    Y0FromTPC( verticalAlignment[eVTPC2vsMTPCR], y0shifts[eVTPC2vsMTPCL], y0shifts[eMTPCRvsVTPC2_crosscheck],true );
    y0shifts[eMTPCRvsVTPC2_crosscheck].tpcName = "MTPCR_crosscheck";
    //GTPC y shift crosscheck from VTPC2
    Y0FromTPC( verticalAlignment[eVTPC2vsMTPCL], y0shifts[eVTPC2vsMTPCL], y0shifts[eGTPCvsVTPC2]);
    y0shifts[eGTPCvsVTPC2].tpcName = "GTPC_crosscheck";

    /*
    PrintY0Shifts( y0shifts[eMTPCLvsTOFL] );
    PrintY0Shifts( y0shifts[eMTPCRvsTOFR] );
    PrintY0Shifts( y0shifts[eVTPC2vsMTPCL] );
    PrintY0Shifts( y0shifts[eVTPC1vsVTPC2] );
    PrintY0Shifts( y0shifts[eVTPC1vsGTPC] );
    PrintY0Shifts( y0shifts[eMTPCRvsVTPC2_crosscheck] );
    PrintY0Shifts( y0shifts[eGTPCvsVTPC2] );
    */
//     int timeOffset = 0;
//     if(y0shifts[eVTPC1vsVTPC2].sections.size()!=0 )
//     {
//       timeOffset = y0shifts[eVTPC1vsVTPC2].sections[0].starttime;
//     }

    for(std::map<eCommonPlanes, y0_corrections >::iterator it = y0shifts.begin();
        it!= y0shifts.end(); it++)
    {

      PlotY0Shifts( it->second, timeOffset);

      if(!prettyPlot) {
        TCanvas * temp = new TCanvas(("c_y0shift_"+ it->second.tpcName).c_str(),(it->second.tpcName + " vertical shift").c_str(), 1 );
        it->second.g_y0_correction->Draw("ap");
       	if(it->second.smoothFunc) { it->second.smoothFunc->Draw("l same");  }
	temp->Update();
        temp->Print(("c_y0shift_"+ it->second.tpcName + ".eps").c_str());
      }
      std::cout << "Plot for " << it->second.tpcName << " is ready." << std::endl;
      if(it->second.linearFits.size() !=0 ) {
	WriteY0Alignment(it->second);
      }
      else {
	std::cout << "No valid y0 sections for " << it->second.tpcName << ", the calibration file will not be written." << std::endl;
      }
    }
    WriteAbsY0AlignmentXML("TPCAlignmentGeometry.xml", y0shifts);


    if(prettyPlot) {
      PrettyPlot(y0shifts[eMTPCLvsTOFL].g_y0_correction,
		 y0shifts[eMTPCRvsTOFR].g_y0_correction,
		 y0shifts[eVTPC2vsMTPCL].g_y0_correction,
		 y0shifts[eVTPC1vsVTPC2].g_y0_correction,
		 y0shifts[eVTPC1vsGTPC].g_y0_correction,
		 y0shifts[eMTPCRvsVTPC2_crosscheck].g_y0_correction,
		 y0shifts[eGTPCvsVTPC2].g_y0_correction,
		 "vertical shift",
		 "y0shift",
		 y0shifts[eMTPCLvsTOFL].smoothFunc,
		 y0shifts[eMTPCRvsTOFR].smoothFunc,
		 y0shifts[eVTPC2vsMTPCL].smoothFunc,
		 y0shifts[eVTPC1vsVTPC2].smoothFunc,
		 y0shifts[eVTPC1vsGTPC].smoothFunc,
		 y0shifts[eMTPCRvsVTPC2_crosscheck].smoothFunc,
		 y0shifts[eGTPCvsVTPC2].smoothFunc
		 );
    }
  }


  for(std::map<eCommonPlanes, alignmentparams>::iterator it = verticalAlignment.begin();
      it!= verticalAlignment.end(); it++)
  {
    std::string uniqueName_ = it->second.s_TPCvsTPC+"_"+it->second.s_slave+"_"+it->second.s_master;
    std::string uniqueName = it->second.s_TPCvsTPC+" "+it->second.s_slave+" "+it->second.s_master;
    PlotRelativeAlignment(it->second, timeOffset );


    if(!prettyPlot) {
      TCanvas * temp1 = new TCanvas(("c_relative_alignment_"+ uniqueName_ + "slope").c_str(),
                                  (uniqueName + " slope").c_str(), 1 );
      it->second.g_slope->Draw("ap");

      TCanvas * temp2 = new TCanvas(("c_relative_alignment_"+ uniqueName_ + "offset").c_str(),
                                  (uniqueName + " offset").c_str(), 1 );
      it->second.g_offset->Draw("ap");
    }
  }

  if(prettyPlot) {
    PrettyPlot(verticalAlignment[eMTPCLvsTOFL].g_slope,
              verticalAlignment[eMTPCRvsTOFR].g_slope,
              verticalAlignment[eVTPC2vsMTPCL].g_slope,
              verticalAlignment[eVTPC1vsVTPC2].g_slope,
              verticalAlignment[eVTPC1vsGTPC].g_slope,
              verticalAlignment[eVTPC2vsMTPCR].g_slope,
              verticalAlignment[eGTPCvsVTPC2].g_slope,
              "relative alignment slope",
              "relativeAlignment_slope");
   PrettyPlot(verticalAlignment[eMTPCLvsTOFL].g_offset,
              verticalAlignment[eMTPCRvsTOFR].g_offset,
              verticalAlignment[eVTPC2vsMTPCL].g_offset,
              verticalAlignment[eVTPC1vsVTPC2].g_offset,
              verticalAlignment[eVTPC1vsGTPC].g_offset,
              verticalAlignment[eVTPC2vsMTPCR].g_offset,
              verticalAlignment[eGTPCvsVTPC2].g_offset,
              "relative alignment offset",
              "relativeAlignment_offset");
  }
}


int main( int argc, char **argv )
{

  startRange = 0;
  stopRange = 0;

  std::string fileName = "";
  std::string programName = "alignmentVdriftT0";
  bool paramHelp = false;
  bool driftvelocityCampaign = false;
  bool y0Campaign = false;
  bool useRecVdrift = false;
  bool batchMode = false;

  if ( argc<=1 ) {
    paramHelp = true;
    std::cerr << "\nError: Insufficient number of arguments." << std::endl;
  }

  for(int i=1; i!= argc; i++) {
    if(std::string(argv[i]) != "-i" &&
       std::string(argv[i]) != "-v" &&
       std::string(argv[i]) != "-y" &&
       std::string(argv[i]) != "-r" &&
       std::string(argv[i]) != "-a" &&
       std::string(argv[i]) != "-b" &&
       std::string(argv[i]) != "-h")
    {
      paramHelp = true;
      std::cerr << "\nError: Unknown parameter " <<  argv[i] << std::endl;
      break;
    }

    if(std::string(argv[i]) == "-i" && i!=argc-1) {
      fileName = argv[++i];
      continue;
    }
    if(std::string(argv[i]) == "-r") {
      if(i+2>=argc) {
        std::cerr << "Error: Insufficient number of argument to set range." << std::endl;
        paramHelp = true;
        break;
      }
      startRange = atof(argv[++i]);
      stopRange  = atof(argv[++i]);

      if(startRange == stopRange) {
        std::cerr << "Error: StartRange = StopRange" << std::endl;
        paramHelp = true;
        break;
      }
      if(startRange > stopRange) {
        std::cerr << "Error: StartRange > StopRange" << std::endl;
        paramHelp = true;
        break;
      }

      continue;
    }

    if(std::string(argv[i]) == "-b") {
      batchMode = true;
      continue;
    }

    if(std::string(argv[i]) == "-v") {
      driftvelocityCampaign = true;
      continue;
    }

    if(std::string(argv[i]) == "-y") {
      y0Campaign = true;
      continue;
    }

    if(std::string(argv[i]) == "-h") {
      paramHelp = true;
      continue;
    }

    if(std::string(argv[i]) == "-a") {
      useRecVdrift = true;
      continue;
    }

    if(std::string(argv[i]) == "-i" &&
       i == argc-1) {
       paramHelp = true;
       std::cerr << "\nError: Filename was not specified." << std::endl;
       break;
    }
  }

  if(fileName == ""){
    std::cerr << "Error! Please specify input file name" << std::endl;
    paramHelp = true;
  }

  if(!driftvelocityCampaign && !y0Campaign)
  {
    std::cerr << "Error! No campaign were specified" << std::endl;
    paramHelp = true;
  }

  if(paramHelp)
  {
    std::cerr << std::endl
              << "The program " << programName << " is intended to calculate drift velocities and aligment of" << std::endl
              << "the various TPCs, by using the TOF walls as a reference by minimizing the shifts of track" << std::endl
              << "segments matchings between the TPCs and TOF hits." << std::endl << std::endl
              << "Usage:" << std::endl
              << "    " << programName << " <options>" << std::endl << std::endl
              << "Options:" << std::endl
              << "     -i inputFile        single input file name" << std::endl
              << "     -v                  selects drift velocity campaign, assuming good alignment" << std::endl
              << "     -y                  selects vertical alignment campaign," << std::endl
	      << "                         assuming good drift velocities" << std::endl
              << "     -r start stop       plots a specific range (hours)" << std::endl
	      << "     -a                  use drift time instead of cluster y " << std::endl
	      << "                         position and calibrate absolute drift" << std::endl
	      << "                         velocity (experimental, mainTPC only)" << std::endl
	      << "     -b                  batch mode" << std::endl
              << "     -h                  this help message" << std::endl << std::endl
              << "Example:" << std::endl
              << "    " << programName << " -i timeTagged.root -v -y" << std::endl
              << std::endl;
    return -1;
  }

  int dummyArgc = 0;
  char **dummyArgv = 0;
  TApplication calibrator("AlignmentT0Vdrift", &dummyArgc, dummyArgv);
  gROOT->SetStyle("Plain");
  //gStyle->SetOptTitle(kFALSE);
  // gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1);
  gROOT->SetBatch(batchMode); //false by default

/*
  TFile * file = TFile::Open("timetagged.root", "ro");
  alignmentparams mtpcl_tofl;
  mtpcl_tofl.s_slave = slave;
  mtpcl_tofl.s_master = master;

  TTree * dataTree = (TTree *)file->Get(treename.c_str());

  fitVdrift(dataTree, mtpcl_tofl, true, treename);

  TCanvas * c1 = new TCanvas("c1", "tofY", 1);
  mtpcl_tofl.h_master->Draw();

  TCanvas * c2 = new TCanvas("c2", "extrapolated", 1);
  mtpcl_tofl.h_slave->Draw();

  TCanvas * c3 = new TCanvas("c3", "residual", 1);
  mtpcl_tofl.h_residual->Draw();

  TCanvas * c4 = new TCanvas("c4", "fit", 1);
  mtpcl_tofl.h_residual_vs_slave->Draw();

  TCanvas * c5 = new TCanvas("c5", "slope", 1);
  mtpcl_tofl.g_slope->Draw("ap");

  TCanvas * c6 = new TCanvas("c6", "offset", 1);
  mtpcl_tofl.g_offset->Draw("ap");



  TF1 * func = new TF1("m_func","[0]+[1]*x", -70, 70 );
  func->SetParameter(0, mtpcl_tofl.sections[0].offset);
  func->SetParameter(1, mtpcl_tofl.sections[0].slope);
  func->Draw("same");

  std::cout << " f(x) = a+b*x;   a = " << mtpcl_tofl.sections[0].offset << " +- " << mtpcl_tofl.sections[0].cov_offset << " b = " << mtpcl_tofl.sections[0].slope << " +-" << mtpcl_tofl.sections[0].cov_slope << std::endl;


  c1->Print("c1.eps", "eps");
  c2->Print("c2.eps", "eps");
  c3->Print("c3.eps", "eps");
  c4->Print("c4.eps", "eps");
  c5->Print("c5.eps", "eps");
  c6->Print("c6.eps", "eps");
*/

  if(useRecVdrift) {
    AbsDriftVelocityCampaign(fileName);
    //DriftVelocityY0CalibrationCampaign(fileName, driftvelocityCampaign, y0Campaign, useRecVdrift);
  }
  else {
    DriftVelocityY0CalibrationCampaign(fileName, driftvelocityCampaign, y0Campaign, useRecVdrift);
  }
  if(!batchMode) {
    calibrator.Run();
  }
  return 0;
}
