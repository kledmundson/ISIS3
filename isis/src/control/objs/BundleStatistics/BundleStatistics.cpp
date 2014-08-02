#include "BundleStatistics.h"

#include <QDebug>
#include <QDataStream>
#include <QObject>
#include <QString>
#include <QUuid>
#include <QXmlStreamWriter>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Camera.h"
#include "ControlMeasure.h"
#include "ControlNet.h"
#include "ControlPoint.h"
#include "CorrelationMatrix.h"
#include "MaximumLikelihoodWFunctions.h"
#include "PvlKeyword.h"
#include "PvlObject.h"
#include "SerialNumberList.h"
#include "StatCumProbDistDynCalc.h"
#include "Statistics.h"
#include "XmlStackedHandlerReader.h"

using namespace boost::numeric::ublas;
namespace Isis {

  BundleStatistics::BundleStatistics(QObject *parent) : QObject(parent) {

    m_id = NULL;
    m_id = new QUuid(QUuid::createUuid());

    m_correlationMatrix = NULL;
    m_correlationMatrix = new CorrelationMatrix();

    m_numberFixedPoints = 0; // set in BA constructor->init->fillPointIndexMap
    m_numberIgnoredPoints = 0; // set in BA constructor->init->fillPointIndexMap


    // set in BundleAdjust init()
    m_numberHeldImages = 0;

    // members set while computing bundle stats
    m_rmsImageSampleResiduals.clear();
    m_rmsImageLineResiduals.clear();
    m_rmsImageResiduals.clear();
    m_rmsImageXSigmas.clear();
    m_rmsImageYSigmas.clear();
    m_rmsImageZSigmas.clear();
    m_rmsImageRASigmas.clear();
    m_rmsImageDECSigmas.clear();
    m_rmsImageTWISTSigmas.clear();

    m_minSigmaLatitude = 1.0e+12;
    m_minSigmaLatitudePointId = "";

    m_maxSigmaLatitude = 0.0;
    m_maxSigmaLatitudePointId = "";

    m_minSigmaLongitude = 1.0e+12;
    m_minSigmaLongitudePointId = "";

    m_maxSigmaLongitude = 0.0;
    m_maxSigmaLongitudePointId = "";

    m_minSigmaRadius = 1.0e+12;
    m_minSigmaRadiusPointId = "";

    m_maxSigmaRadius = 0.0;
    m_maxSigmaRadiusPointId = "";
    
    m_rmsSigmaLat = 0.0;
    m_rmsSigmaLon = 0.0;
    m_rmsSigmaRad = 0.0;


    // set by compute residuals
    m_rms_rx = 0.0;
    m_rms_ry = 0.0;
    m_rms_rxy = 0.0;

    // set by compute rejection limit
    m_rejectionLimit = 0.0;
    
    // set by flag outliers    
    m_numberRejectedObservations = 0;

    // set by formNormalEquations_CHOLMOD, formNormalEquations_SPECIALK, or solve
    m_numberObservations = 0;
    m_numberImageParameters = 0; // ??? this is the same value an m_nRank
                                 // = observations * numImagePartials

// ??? unused variable ???    m_numberHeldPoints = 0;

// set by formNormalEquations_CHOLMOD, formNormalEquations_SPECIALK, or
// setParameterWeights (i.e. solve)
    m_numberConstrainedPointParameters = 0;
    m_numberConstrainedImageParameters = 0;

    // set by initialize, formNormalEquations_CHOLMOD, formNormalEquations_SPECIALK, or solve
    m_numberUnknownParameters = 0;
    // solve and solve cholesky
    m_degreesOfFreedom = -1;
    m_sigma0 = 0.0;
    m_elapsedTime = 0.0;
    m_elapsedTimeErrorProp = 0.0;
    m_converged = false; // or initialze method

    // SetUp method initializes remaining maximum likelihood estimation parameters:
    // m_wFunc, m_cumPro, m_maxLikelihoodIndex, m_maxLikelihoodQuan,
    // m_maxLikelihoodMedianR2Residuals
    m_cumPro = NULL;
    m_cumProRes = NULL;
    m_numberMaximumLikelihoodModels = 0;
    m_maximumLikelihoodIndex = 0;
    m_maximumLikelihoodMedianR2Residuals = 0.0;
    m_wFunc[0] = m_wFunc[1] = m_wFunc[2] = NULL;
    m_maximumLikelihoodQuan[0] = m_maximumLikelihoodQuan[1] = m_maximumLikelihoodQuan[2] = 0.5; // better init value ???


    // residual prob distribution is calculated even if there is no maximum likelihood estimation
    // set up the solver to have a node at every percent of the distribution
    m_cumProRes = new StatCumProbDistDynCalc;
    initializeResidualsProbabilityDistribution(101);

  }



  /**
   * Construct this BundleSettings object from XML.
   *
   * @param bundleSettingsFolder Where this settings XML resides - /work/.../projectRoot/images/import1
   * @param xmlReader An XML reader that's up to an <bundleSettings/> tag.
   * @param parent The Qt-relationship parent
   */
  BundleStatistics::BundleStatistics(Project *project, XmlStackedHandlerReader *xmlReader,
                                 QObject *parent) : QObject(parent) {
    m_id = NULL;
    // ??? initializations ???
    xmlReader->pushContentHandler(new XmlHandler(this, project));
  }



  BundleStatistics::BundleStatistics(const BundleStatistics &other)
      : m_id(new QUuid(other.m_id->toString())),
        m_correlationMatrix(new CorrelationMatrix(*other.m_correlationMatrix)),
        m_numberFixedPoints(other.m_numberFixedPoints),
        m_numberIgnoredPoints(other.m_numberIgnoredPoints),
        m_numberHeldImages(other.m_numberHeldImages),
        m_rms_rx(other.m_rms_rx),
        m_rms_ry(other.m_rms_ry),
        m_rms_rxy(other.m_rms_rxy),
        m_rejectionLimit(other.m_rejectionLimit),
        m_numberRejectedObservations(other.m_numberRejectedObservations),
        m_numberObservations(other.m_numberObservations),
        m_numberImageParameters(other.m_numberImageParameters),
        m_numberConstrainedPointParameters(other.m_numberConstrainedPointParameters),
        m_numberConstrainedImageParameters(other.m_numberConstrainedImageParameters),
        m_numberUnknownParameters(other.m_numberUnknownParameters),
        m_degreesOfFreedom(other.m_degreesOfFreedom),
        m_sigma0(other.m_sigma0),
        m_elapsedTime(other.m_elapsedTime),
        m_elapsedTimeErrorProp(other.m_elapsedTimeErrorProp),
        m_converged(other.m_converged),
        m_rmsImageSampleResiduals(other.m_rmsImageSampleResiduals),
        m_rmsImageLineResiduals(other.m_rmsImageLineResiduals),
        m_rmsImageResiduals(other.m_rmsImageResiduals),
        m_rmsImageXSigmas(other.m_rmsImageXSigmas),    
        m_rmsImageYSigmas(other.m_rmsImageYSigmas),    
        m_rmsImageZSigmas(other.m_rmsImageZSigmas),    
        m_rmsImageRASigmas(other.m_rmsImageRASigmas),   
        m_rmsImageDECSigmas(other.m_rmsImageDECSigmas),  
        m_rmsImageTWISTSigmas(other.m_rmsImageTWISTSigmas),
        m_minSigmaLatitude(other.m_minSigmaLatitude),
        m_minSigmaLatitudePointId(other.m_minSigmaLatitudePointId),
        m_maxSigmaLatitude(other.m_maxSigmaLatitude),
        m_maxSigmaLatitudePointId(other.m_maxSigmaLatitudePointId),
        m_minSigmaLongitude(other.m_minSigmaLongitude),
        m_minSigmaLongitudePointId(other.m_minSigmaLongitudePointId),
        m_maxSigmaLongitude(other.m_maxSigmaLongitude),
        m_maxSigmaLongitudePointId(other.m_maxSigmaLongitudePointId),
        m_minSigmaRadius(other.m_minSigmaRadius),
        m_minSigmaRadiusPointId(other.m_minSigmaRadiusPointId),
        m_maxSigmaRadius(other.m_maxSigmaRadius),
        m_maxSigmaRadiusPointId(other.m_maxSigmaRadiusPointId),
        m_rmsSigmaLat(other.m_rmsSigmaLat),
        m_rmsSigmaLon(other.m_rmsSigmaLon),
        m_rmsSigmaRad(other.m_rmsSigmaRad),
        m_numberMaximumLikelihoodModels(other.m_numberMaximumLikelihoodModels),
        m_maximumLikelihoodIndex(other.m_maximumLikelihoodIndex),
        m_cumPro(new StatCumProbDistDynCalc(*other.m_cumPro)),
        m_cumProRes(new StatCumProbDistDynCalc(*other.m_cumProRes)), 
        m_maximumLikelihoodMedianR2Residuals(other.m_maximumLikelihoodMedianR2Residuals) {

    for (int i = 0; i < 3; i++) {
      m_wFunc[i] = new MaximumLikelihoodWFunctions(*other.m_wFunc[i]);
      m_maximumLikelihoodQuan[i] = other.m_maximumLikelihoodQuan[i];
    }

  }



  BundleStatistics::~BundleStatistics() {
    
    delete m_id;
    m_id = NULL;

    delete m_correlationMatrix;
    m_correlationMatrix = NULL;

    delete m_cumPro;
    m_cumPro = NULL;

    delete m_cumProRes;
    m_cumProRes = NULL;

    for (int i = 0; i < m_numberMaximumLikelihoodModels; i++) {
      delete m_wFunc[i];
      m_wFunc[i] = NULL;
    }

  }

  

  BundleStatistics &BundleStatistics::operator=(const BundleStatistics &other) {

    if (&other != this) {
      delete m_id;
      m_id = NULL;
      m_id = new QUuid(m_id->toString());

      delete m_correlationMatrix;
      m_correlationMatrix = NULL;
      m_correlationMatrix = new CorrelationMatrix(*other.m_correlationMatrix);

      m_numberFixedPoints = other.m_numberFixedPoints;
      m_numberIgnoredPoints = other.m_numberIgnoredPoints;
      m_numberHeldImages = other.m_numberHeldImages;
      m_rms_rx = other.m_rms_rx;
      m_rms_ry = other.m_rms_ry;
      m_rms_rxy = other.m_rms_rxy;
      m_rejectionLimit = other.m_rejectionLimit;
      m_numberRejectedObservations = other.m_numberRejectedObservations;
      m_numberObservations = other.m_numberObservations;
      m_numberImageParameters = other.m_numberImageParameters;
      m_numberConstrainedPointParameters = other.m_numberConstrainedPointParameters;
      m_numberConstrainedImageParameters = other.m_numberConstrainedImageParameters;
      m_numberUnknownParameters = other.m_numberUnknownParameters;
      m_degreesOfFreedom = other.m_degreesOfFreedom;
      m_sigma0 = other.m_sigma0;
      m_elapsedTime = other.m_elapsedTime;
      m_elapsedTimeErrorProp = other.m_elapsedTimeErrorProp;
      m_converged = other.m_converged;
      m_rmsImageSampleResiduals = other.m_rmsImageSampleResiduals;
      m_rmsImageLineResiduals = other.m_rmsImageLineResiduals;
      m_rmsImageResiduals = other.m_rmsImageResiduals;
      m_rmsImageXSigmas = other.m_rmsImageXSigmas;
      m_rmsImageYSigmas = other.m_rmsImageYSigmas;
      m_rmsImageZSigmas = other.m_rmsImageZSigmas;
      m_rmsImageRASigmas = other.m_rmsImageRASigmas;
      m_rmsImageDECSigmas = other.m_rmsImageDECSigmas;
      m_rmsImageTWISTSigmas = other.m_rmsImageTWISTSigmas;
      m_minSigmaLatitude = other.m_minSigmaLatitude;
      m_minSigmaLatitudePointId = other.m_minSigmaLatitudePointId;
      m_maxSigmaLatitude = other.m_maxSigmaLatitude;
      m_maxSigmaLatitudePointId = other.m_maxSigmaLatitudePointId;
      m_minSigmaLongitude = other.m_minSigmaLongitude;
      m_minSigmaLongitudePointId = other.m_minSigmaLongitudePointId;
      m_maxSigmaLongitude = other.m_maxSigmaLongitude;
      m_maxSigmaLongitudePointId = other.m_maxSigmaLongitudePointId;
      m_minSigmaRadius = other.m_minSigmaRadius;
      m_minSigmaRadiusPointId = other.m_minSigmaRadiusPointId;
      m_maxSigmaRadius = other.m_maxSigmaRadius;
      m_maxSigmaRadiusPointId = other.m_maxSigmaRadiusPointId;
      m_rmsSigmaLat = other.m_rmsSigmaLat;
      m_rmsSigmaLon = other.m_rmsSigmaLon;
      m_rmsSigmaRad = other.m_rmsSigmaRad;
      m_numberMaximumLikelihoodModels = other.m_numberMaximumLikelihoodModels;

      for (int i = 0; i < m_numberMaximumLikelihoodModels; i++) {

        delete m_wFunc[i];
        m_wFunc[i] = new MaximumLikelihoodWFunctions(*other.m_wFunc[i]);

        m_maximumLikelihoodQuan[i] = other.m_maximumLikelihoodQuan[i];
      }

      m_maximumLikelihoodIndex = other.m_maximumLikelihoodIndex;

      delete m_cumPro;
      m_cumPro = new StatCumProbDistDynCalc(*other.m_cumPro);

      delete m_cumProRes;
      m_cumProRes = new StatCumProbDistDynCalc(*other.m_cumProRes);

      m_maximumLikelihoodMedianR2Residuals = other.m_maximumLikelihoodMedianR2Residuals;
    }
    return *this;
  }




  /**
   * Compute Bundle statistics. 
   *  
   * Sets: 
   * m_rmsImageSampleResiduals 
   * m_rmsImageLineResiduals 
   * m_rmsImageResiduals 
   * 
   * m_rmsImageXSigmas     UNUSED???
   * m_rmsImageYSigmas     UNUSED???
   * m_rmsImageZSigmas     UNUSED???
   * m_rmsImageRASigmas    UNUSED???
   * m_rmsImageDECSigmas   UNUSED???
   * m_rmsImageTWISTSigmas UNUSED???
   *  
   * m_maxSigmaLatitude 
   * m_maxSigmaLatitudePointId 
   * m_maxSigmaLongitude 
   * m_maxSigmaLongitudePointId 
   * m_maxSigmaRadius 
   * m_maxSigmaRadiusPointId
   *                       
   * m_minSigmaLatitude   
   * m_minSigmaLatitudePointId  
   * m_minSigmaLongitude  
   * m_minSigmaLongitudePointId 
   * m_minSigmaRadius     
   * m_minSigmaRadiusPointId    
   *  
   * m_rmsSigmaLat
   * m_rmsSigmaLon
   * m_rmsSigmaRad
   *  
   *  
   *  Needed:
   * SerialNumberList *snList, ControlNet *cnet, bool errorPropagation,
   * bool solveRadius 
   *
   *
   */
  bool BundleStatistics::computeBundleStatistics(SerialNumberList *snList, 
                                                 ControlNet *cnet, 
                                                 bool errorPropagation,
                                                 bool solveRadius) {
    int    numImages      = snList->Size();
    int    numMeasures    = 0;
    int    imageIndex     = 0;
    double sampleResidual = 0;
    double lineResidual   = 0;

    m_rmsImageSampleResiduals.resize(numImages);
    m_rmsImageLineResiduals.resize(numImages);
    m_rmsImageResiduals.resize(numImages);

    // load image coordinate residuals into statistics vectors
    int numObs = 0;
    int numObjectPoints = cnet->GetNumPoints();

    for (int i = 0; i < numObjectPoints; i++) {

      const ControlPoint *point = cnet->GetPoint(i);
      if (point->IsIgnored()) {
        continue;
      }

      if (point->IsRejected()) {
        continue;
      }

      numMeasures = point->GetNumMeasures();
      for (int j = 0; j < numMeasures; j++) {

        const ControlMeasure *measure = point->GetMeasure(j);
        if (measure->IsIgnored()) {
          continue;
        }

        if (measure->IsRejected()) {
          continue;
        }

        sampleResidual = fabs(measure->GetSampleResidual());
        lineResidual = fabs(measure->GetLineResidual());

        // Determine the image index
        imageIndex = snList->SerialNumberIndex(measure->GetCubeSerialNumber());

        // add residuals to pertinent vector
        m_rmsImageSampleResiduals[imageIndex].AddData(sampleResidual);
        m_rmsImageLineResiduals[imageIndex].AddData(lineResidual);
        m_rmsImageResiduals[imageIndex].AddData(lineResidual);
        m_rmsImageResiduals[imageIndex].AddData(sampleResidual);

        numObs++;
      }
    }

    if (errorPropagation) {
      m_rmsImageXSigmas.resize(numImages);
      m_rmsImageYSigmas.resize(numImages);
      m_rmsImageZSigmas.resize(numImages);
      m_rmsImageRASigmas.resize(numImages);
      m_rmsImageDECSigmas.resize(numImages);
      m_rmsImageTWISTSigmas.resize(numImages);

      // compute stats for point sigmas
      Statistics sigmaLatitude;
      Statistics sigmaLongitude;
      Statistics sigmaRadius;

      double     dSigmaLat, dSigmaLong, dSigmaRadius;

      int        nPoints        = cnet->GetNumPoints();
      for (int i = 0; i < nPoints; i++) {

        const ControlPoint *point = cnet->GetPoint(i);
        if (point->IsIgnored()) {
          continue;
        }

        dSigmaLat    = point->GetAdjustedSurfacePoint().GetLatSigmaDistance().meters();
        dSigmaLong   = point->GetAdjustedSurfacePoint().GetLonSigmaDistance().meters();
        dSigmaRadius = point->GetAdjustedSurfacePoint().GetLocalRadiusSigma().meters();

        sigmaLatitude.AddData(dSigmaLat);
        sigmaLongitude.AddData(dSigmaLong);
        sigmaRadius.AddData(dSigmaRadius);

        if (i > 0) {
          if (dSigmaLat > m_maxSigmaLatitude) {
            m_maxSigmaLatitude = dSigmaLat;
            m_maxSigmaLatitudePointId = point->GetId();
          }
          if (dSigmaLong > m_maxSigmaLongitude) {
            m_maxSigmaLongitude = dSigmaLong;
            m_maxSigmaLongitudePointId = point->GetId();
          }
          if (solveRadius) {
            if (dSigmaRadius > m_maxSigmaRadius) {
              m_maxSigmaRadius = dSigmaRadius;
              m_maxSigmaRadiusPointId = point->GetId();
            }
          }
          if (dSigmaLat < m_minSigmaLatitude) {
            m_minSigmaLatitude = dSigmaLat;
            m_minSigmaLatitudePointId = point->GetId();
          }
          if (dSigmaLong < m_minSigmaLongitude) {
            m_minSigmaLongitude = dSigmaLong;
            m_minSigmaLongitudePointId = point->GetId();
          }
          if (solveRadius) {
            if (dSigmaRadius < m_minSigmaRadius) {
              m_minSigmaRadius = dSigmaRadius;
              m_minSigmaRadiusPointId = point->GetId();
            }
          }
        }
        else {
          m_maxSigmaLatitude = dSigmaLat;
          m_maxSigmaLongitude = dSigmaLong;

          m_minSigmaLatitude = dSigmaLat;
          m_minSigmaLongitude = dSigmaLong;

          if (solveRadius) {
            m_maxSigmaRadius = dSigmaRadius;
            m_minSigmaRadius = dSigmaRadius;
          }
        }
      }
      m_rmsSigmaLat = sigmaLatitude.Rms();
      m_rmsSigmaLon = sigmaLongitude.Rms();
      m_rmsSigmaRad = sigmaRadius.Rms();
    }

    return true;
  }



  /** 
  * This method steps up the maximum likelihood estimation solution.  Up to three successive
  * solutions models are available.
  */
  void BundleStatistics::maximumLikelihoodSetUp(
      QList< QPair< MaximumLikelihoodWFunctions::Model, double > > modelsWithQuantiles) {

    // reinitialize variables if this setup has already been called
    if (m_numberMaximumLikelihoodModels > 0) {

      m_numberMaximumLikelihoodModels = 0;
      m_maximumLikelihoodIndex = 0;
      m_maximumLikelihoodMedianR2Residuals = 0.0;
      m_maximumLikelihoodQuan[0] = m_maximumLikelihoodQuan[1] = m_maximumLikelihoodQuan[2] = 0.5; // better init value ???

      delete m_cumPro;
      m_cumPro = NULL;

      delete m_cumProRes;
      m_cumProRes = NULL;


      for (int i = 0;i < m_numberMaximumLikelihoodModels;i++) {
        delete m_wFunc[i];
        m_wFunc[i] = NULL;
      }

    }

    // initialize to NULL
    // JWB - potentially, we could generalize for any size QList
    m_numberMaximumLikelihoodModels = modelsWithQuantiles.size();

    // residual prob distribution is calculated even if there is no maximum likelihood estimation
    // set up the solver to have a node at every percent of the distribution
    m_cumProRes = new StatCumProbDistDynCalc;
    initializeResidualsProbabilityDistribution(101);

    // these values are only set if there will be maximum likelihood estimation...
    if (modelsWithQuantiles.size() != 0) {
      // MaximumLikeliHood Estimation is being used

      // set up the cumulative prob solver to have a node at every percent of the distribution
      m_cumPro = new StatCumProbDistDynCalc;
      initializeProbabilityDistribution(101);

      // set up the w functions
      for (int i = 0; i < m_numberMaximumLikelihoodModels; i++) {

        m_wFunc[i] = new MaximumLikelihoodWFunctions;
        m_wFunc[i]->setModel(modelsWithQuantiles[i].first);

        m_maximumLikelihoodQuan[i] = modelsWithQuantiles[i].second;

      }
    }

    //maximum likelihood estimation tiered solutions requiring multiple convergeances are support,
    // this index keeps track of which tier the solution is in
    m_maximumLikelihoodIndex = 0;
  }



  void BundleStatistics::printMaximumLikelihoodTierInformation() {
    printf("Maximum Likelihood Tier: %d\n", m_maximumLikelihoodIndex);
    if (m_numberMaximumLikelihoodModels > m_maximumLikelihoodIndex) {
      // if maximum likelihood estimation is being used
      // at the end of every iteration
      // reset the tweaking contant to the desired quantile of the |residual| distribution
      m_wFunc[m_maximumLikelihoodIndex]->setTweakingConstant(m_cumPro->value(m_maximumLikelihoodQuan[m_maximumLikelihoodIndex]));
      //  print meadians of residuals
      m_maximumLikelihoodMedianR2Residuals = m_cumPro->value(0.5);
      printf("Median of R^2 residuals:  %lf\n", m_maximumLikelihoodMedianR2Residuals);
      //restart the dynamic calculation of the cumulative probility distribution of |R^2 residuals| --so it will be up to date for the next iteration
      initializeProbabilityDistribution(101);
    }
  }



  void BundleStatistics::initializeProbabilityDistribution(unsigned int nodes) {
    m_cumPro->initialize(nodes);
  }



  void BundleStatistics::initializeResidualsProbabilityDistribution(unsigned int nodes) {
    m_cumProRes->initialize(nodes);
  }



  void BundleStatistics::addProbabilityDistributionObservation(double observationValue) {
    m_cumPro->addObs(observationValue);
  }



  void BundleStatistics::addResidualsProbabilityDistributionObservation(double observationValue) {
    m_cumProRes->addObs(observationValue);
  }



  void BundleStatistics::incrementMaximumLikelihoodModelIndex() {
    m_maximumLikelihoodIndex++;
  }



  void BundleStatistics::incrementFixedPoints() {
    m_numberFixedPoints++;
  }



  int BundleStatistics::numberFixedPoints() {
    return m_numberFixedPoints;
  }



  void BundleStatistics::incrementHeldImages() {
    m_numberHeldImages++;
  }



  int BundleStatistics::numberHeldImages() {
    return m_numberHeldImages;
  }



  void BundleStatistics::incrementIgnoredPoints() {
    m_numberIgnoredPoints++;
  }


  void BundleStatistics::setRmsXYResiduals(double rx, double ry, double rxy) {
    m_rms_rx = rx;
    m_rms_ry = ry;
    m_rms_rxy = rxy;
  }



  void BundleStatistics::setRmsRx(double rx) {
    m_rms_rx = rx;
  }



  void BundleStatistics::setRmsRy(double ry) {
    m_rms_ry = ry;
  }



  void BundleStatistics::setRmsRxy(double rxy) {
    m_rms_rxy = rxy;
  }



  void BundleStatistics::setRejectionLimit(double rejectionLimit) {
    m_rejectionLimit = rejectionLimit;
  }



  void BundleStatistics::setNumberRejectedObservations(int numberRejectedObservations) {
    m_numberRejectedObservations = numberRejectedObservations;
  }



  void BundleStatistics::setNumberObservations(int numberObservations) {
    m_numberObservations = numberObservations;
  }



  void BundleStatistics::setNumberImageParameters(int numberParameters) {
    m_numberImageParameters = numberParameters;
  }



  void BundleStatistics::resetNumberConstrainedPointParameters() {
    m_numberConstrainedPointParameters = 0;
  }



  void BundleStatistics::incrementNumberConstrainedPointParameters(int incrementAmount) {
    m_numberConstrainedPointParameters += incrementAmount;
  }



  void BundleStatistics::resetNumberConstrainedImageParameters() {
    m_numberConstrainedImageParameters = 0;
  }



  void BundleStatistics::incrementNumberConstrainedImageParameters(int incrementAmount) {
    m_numberConstrainedImageParameters += incrementAmount;
  }



  void BundleStatistics::setNumberUnknownParameters(int numberParameters) {
    m_numberUnknownParameters = numberParameters;
  }



  void BundleStatistics::computeDegreesOfFreedom() {
    m_degreesOfFreedom = m_numberObservations
                         + m_numberConstrainedPointParameters
                         + m_numberConstrainedImageParameters
                         - m_numberUnknownParameters;
  }



  void BundleStatistics::computeSigma0(double dvtpv, BundleSettings::ConvergenceCriteria criteria) {
    computeDegreesOfFreedom();

    if (m_degreesOfFreedom > 0) {
      m_sigma0 = dvtpv / m_degreesOfFreedom;
    }
    else if (m_degreesOfFreedom == 0 && criteria == BundleSettings::ParameterCorrections) {
      m_sigma0 = dvtpv;
    }
    else {
      QString msg = "Degrees of Freedom " + toString(m_degreesOfFreedom)
                    + " is invalid (&lt;= 0)!";
      throw IException(IException::Io, msg, _FILEINFO_);
    }

    m_sigma0 = sqrt(m_sigma0);
  }



  void BundleStatistics::setDegreesOfFreedom(double degreesOfFreedom) { // old sparse
    m_degreesOfFreedom = degreesOfFreedom;
  }



  void BundleStatistics::setSigma0(double sigma0) { // old sparse
    m_sigma0 = sigma0;
  }



  void BundleStatistics::setElapsedTime(double time) {
    m_elapsedTime = time;
  }



  void BundleStatistics::setElapsedTimeErrorProp(double time) {
    m_elapsedTimeErrorProp = time;
  }



  void BundleStatistics::setConverged(bool converged) {
    m_converged = converged;
  }



  //************************* Accessors **********************************************************//
  QVector< Statistics > BundleStatistics::rmsImageSampleResiduals() const {
    return m_rmsImageSampleResiduals;
  }



  QVector< Statistics > BundleStatistics::rmsImageLineResiduals() const {
    return m_rmsImageLineResiduals;
  }



  QVector< Statistics > BundleStatistics::rmsImageResiduals() const {
    return m_rmsImageResiduals;
  }



  // currently unused ??? QVector< Statistics > rmsImageXSigmas() const;       
  // currently unused ??? QVector< Statistics > rmsImageYSigmas() const;       
  // currently unused ??? QVector< Statistics > rmsImageZSigmas() const;       
  // currently unused ??? QVector< Statistics > rmsImageRASigmas() const;      
  // currently unused ??? QVector< Statistics > rmsImageDECSigmas() const;     
  // currently unused ??? QVector< Statistics > rmsImageTWISTSigmas() const;   



  double  BundleStatistics::minSigmaLatitude() const {
    return m_minSigmaLatitude;
  }



  QString BundleStatistics::minSigmaLatitudePointId() const {
    return m_minSigmaLatitudePointId;
  }



  double  BundleStatistics::maxSigmaLatitude() const {
    return m_maxSigmaLatitude;
  }



  QString BundleStatistics::maxSigmaLatitudePointId() const {
    return m_maxSigmaLatitudePointId;
  }



  double  BundleStatistics::minSigmaLongitude() const {
    return m_minSigmaLongitude;
  }



  QString BundleStatistics::minSigmaLongitudePointId() const {
    return m_minSigmaLongitudePointId;
  }



  double  BundleStatistics::maxSigmaLongitude() const {
    return m_maxSigmaLongitude;
  }



  QString BundleStatistics::maxSigmaLongitudePointId() const {
    return m_maxSigmaLongitudePointId;
  }



  double  BundleStatistics::minSigmaRadius() const {
    return m_minSigmaRadius;
  }



  QString BundleStatistics::minSigmaRadiusPointId() const {
    return m_minSigmaRadiusPointId;
  }



  double  BundleStatistics::maxSigmaRadius() const {
    return m_maxSigmaRadius;
  }



  QString BundleStatistics::maxSigmaRadiusPointId() const {
    return m_maxSigmaRadiusPointId;
  }



  double BundleStatistics::rmsSigmaLat() const {
    return m_rmsSigmaLat;
  }



  double BundleStatistics::rmsSigmaLon() const {
    return m_rmsSigmaLon;
  }



  double BundleStatistics::rmsSigmaRad() const {
    return m_rmsSigmaRad;
  }



  // currently unused ??? double BundleStatistics::rmsRx() const;
  // currently unused ??? double BundleStatistics::rmsRy() const;
  // currently unused ??? double BundleStatistics::rmsRxy() const; 



  double BundleStatistics::rejectionLimit() const {
    return m_rejectionLimit;
  }



  int BundleStatistics::numberRejectedObservations() const {
    return m_numberRejectedObservations;
  }



  int BundleStatistics::numberObservations() const {
    return m_numberObservations;
  }



  int BundleStatistics::numberImageParameters() const { // ??? this is the same value an m_nRank
    return m_numberImageParameters;
  }



  int BundleStatistics::numberConstrainedPointParameters() const {
    return m_numberConstrainedPointParameters;
  }



  int BundleStatistics::numberConstrainedImageParameters() const {
    return m_numberConstrainedImageParameters;
  }



  int BundleStatistics::numberUnknownParameters() const {
    return m_numberUnknownParameters;
  }



  int BundleStatistics::degreesOfFreedom() const {
    return m_degreesOfFreedom;
  }



  double BundleStatistics::sigma0() const {
    return m_sigma0;
  }



  double BundleStatistics::elapsedTime() const {
    return m_elapsedTime;
  }



  double BundleStatistics::elapsedTimeErrorProp() const {
    return m_elapsedTimeErrorProp;
  }



  bool BundleStatistics::converged() const {
    return m_converged;
  }



  int BundleStatistics::numberMaximumLikelihoodModels() const {
    return m_numberMaximumLikelihoodModels;
  }



  int BundleStatistics::maximumLikelihoodModelIndex() const {
    return m_maximumLikelihoodIndex;
  }



  StatCumProbDistDynCalc BundleStatistics::cumulativeProbabilityDistribution() const {
    return *m_cumPro;
  }



  StatCumProbDistDynCalc BundleStatistics::residualsCumulativeProbabilityDistribution() const {
    return *m_cumProRes;
  }



  double BundleStatistics::maximumLikelihoodMedianR2Residuals() const {
    return m_maximumLikelihoodMedianR2Residuals;
  }



  MaximumLikelihoodWFunctions BundleStatistics::maximumLikelihoodModelWFunc(int modelIndex) const {
    return *m_wFunc[modelIndex];
  }



  double BundleStatistics::maximumLikelihoodModelQuantile(int modelIndex) const {
    return m_maximumLikelihoodQuan[modelIndex];
  }



  PvlObject BundleStatistics::pvlObject(QString name) const {

    PvlObject pvl(name);

    pvl += PvlKeyword("CorrelationMatrix", toString(bool(m_correlationMatrix != NULL)));
    //covariance file name and location
    pvl += PvlKeyword("NumberFixedPoints", toString(m_numberFixedPoints));
    pvl += PvlKeyword("NumberIgnoredPoints", toString(m_numberIgnoredPoints));
    pvl += PvlKeyword("NumberHeldImages", toString(m_numberHeldImages));
    pvl += PvlKeyword("RMSResidualX", toString(m_rms_rx));
    pvl += PvlKeyword("RMSResidualY", toString(m_rms_ry));
    pvl += PvlKeyword("RMSResidualXY", toString(m_rms_rxy));
    pvl += PvlKeyword("RejectionLimit", toString(m_rejectionLimit));
    pvl += PvlKeyword("NumberRejectedObservations", toString(m_numberRejectedObservations));
    pvl += PvlKeyword("NumberObservations", toString(m_numberObservations));
    pvl += PvlKeyword("NumberImageParameters", toString(m_numberImageParameters));
    pvl += PvlKeyword("NumberConstrainedPointParameters", toString(m_numberConstrainedPointParameters));
    pvl += PvlKeyword("NumberConstrainedImageParameters", toString(m_numberConstrainedImageParameters));
    pvl += PvlKeyword("NumberUnknownParameters", toString(m_numberUnknownParameters));
    pvl += PvlKeyword("DegreesOfFreedom", toString(m_degreesOfFreedom));
    pvl += PvlKeyword("Sigma0", toString(m_sigma0));
    pvl += PvlKeyword("ElapsedTime", toString(m_elapsedTime));
    pvl += PvlKeyword("ElapsedTimeErrorProp", toString(m_elapsedTimeErrorProp));
    pvl += PvlKeyword("Converged", toString(m_converged));
    // loop through these ???
    // for (int i = 0; i < m_rmsImageSampleResiduals.size(); i++) {
    //   pvl += PvlKeyword("// QVector< Statistics > rmsImageSampleResiduals", toString(m_));
    // }
    // pvl += PvlKeyword("// QVector< Statistics > rmsImageLineResiduals", toString(m_));
    // pvl += PvlKeyword("// QVector< Statistics > rmsImageResiduals", toString(m_));
    // pvl += PvlKeyword("// QVector< Statistics > rmsImageXSigmas", toString(m_));
    // pvl += PvlKeyword("// QVector< Statistics > rmsImageYSigmas", toString(m_));
    // pvl += PvlKeyword("// QVector< Statistics > rmsImageZSigmas", toString(m_));
    // pvl += PvlKeyword("// QVector< Statistics > rmsImageRASigmas", toString(m_));
    // pvl += PvlKeyword("// QVector< Statistics > rmsImageDECSigmas", toString(m_));
    // pvl += PvlKeyword("// QVector< Statistics > rmsImageTWISTSigmas", toString(m_));
    pvl += PvlKeyword("MinSigmaLatitude", toString(m_minSigmaLatitude));
    pvl += PvlKeyword("MinSigmaLatitudePointId", m_minSigmaLatitudePointId);
    pvl += PvlKeyword("MaxSigmaLatitude", toString(m_maxSigmaLatitude));
    pvl += PvlKeyword("MaxSigmaLatitudePointId", m_maxSigmaLatitudePointId);
    pvl += PvlKeyword("MinSigmaLongitude", toString(m_minSigmaLongitude));
    pvl += PvlKeyword("MinSigmaLongitudePointId", m_minSigmaLongitudePointId);
    pvl += PvlKeyword("MaxSigmaLongitude", toString(m_maxSigmaLongitude));
    pvl += PvlKeyword("MaxSigmaLongitudePointId", m_maxSigmaLongitudePointId);
    pvl += PvlKeyword("MinSigmaRadius", toString(m_minSigmaRadius));
    pvl += PvlKeyword("MinSigmaRadiusPointId", m_minSigmaRadiusPointId);
    pvl += PvlKeyword("MaxSigmaRadius", toString(m_maxSigmaRadius));
    pvl += PvlKeyword("MaxSigmaRadiusPointId", m_maxSigmaRadiusPointId);
    pvl += PvlKeyword("RmsSigmaLat", toString(m_rmsSigmaLat));
    pvl += PvlKeyword("RmsSigmaLon", toString(m_rmsSigmaLon));
    pvl += PvlKeyword("RmsSigmaRad", toString(m_rmsSigmaRad));
    pvl += PvlKeyword("NumberMaximumLikelihoodModels", toString(m_numberMaximumLikelihoodModels));
    if (m_numberMaximumLikelihoodModels > 0) {

      PvlKeyword models("MaximumLikelihoodModels",
                        MaximumLikelihoodWFunctions::modelToString(m_wFunc[0]->model()));
      
      PvlKeyword quantiles("MaximumLikelihoodQuantiles", toString(m_maximumLikelihoodQuan[0]));

      for (int i = 1; i < m_numberMaximumLikelihoodModels; i++) {
        models.addValue(MaximumLikelihoodWFunctions::modelToString(m_wFunc[i]->model()));
        quantiles.addValue(toString(m_maximumLikelihoodQuan[i]));
      }
      pvl += models;
      pvl += quantiles;
      pvl += PvlKeyword("MaximumLikelihoodMedianR2Residuals", 
                          toString(m_maximumLikelihoodMedianR2Residuals));
    }

    return pvl;
  }



  /**
   * Accessor for the Correlation Matrix.
   *
   * @return The correlation matrix.
   */
  CorrelationMatrix BundleStatistics::correlationMatrix() {
    return *m_correlationMatrix;
  }



  /**
   * Set the covariance file name for the matrix used to calculate the correlation matrix.
   *
   * @param name Name of the file used to store the covariance matrix.
   */
   void BundleStatistics::setCorrMatCovFileName(FileName name) {
     m_correlationMatrix->setCovarianceFileName(name);
   }



  /**
   * Set the images and their associated parameters.
   *
   * @param imgsAndParams The qmap with all the images and parameters used for this bundle.
   */
   void BundleStatistics::setCorrMatImgsAndParams(QMap<QString, QStringList> imgsAndParams) {
     m_correlationMatrix->setImagesAndParameters(imgsAndParams);
   }



  void BundleStatistics::save(QXmlStreamWriter &stream, const Project *project, 
                                            FileName newProjectRoot) const {

    stream.writeTextElement("id", m_id->toString());
 
    stream.writeStartElement("instrumentId");
    stream.writeCharacters(m_instrumentId);
    stream.writeEndElement();

    stream.writeStartElement("correlationMatrix");
    // TODO ??? m_correlationMatrix.save(stream, project, newProjectRoot);
    stream.writeEndElement();
    
    stream.writeStartElement("generalStatisticsValues");
    stream.writeStartElement("numberFixedPoints");
    stream.writeCharacters(toString(m_numberFixedPoints));
    stream.writeEndElement();
    stream.writeStartElement("numberIgnoredPoints");
    stream.writeCharacters(toString(m_numberIgnoredPoints));
    stream.writeEndElement();
    stream.writeStartElement("numberHeldImages");
    stream.writeCharacters(toString(m_numberHeldImages));
    stream.writeEndElement();
    stream.writeStartElement("rejectionLimit");
    stream.writeCharacters(toString(m_rejectionLimit));
    stream.writeEndElement();
    stream.writeStartElement("numberRejectedObservations");
    stream.writeCharacters(toString(m_numberRejectedObservations));
    stream.writeEndElement();
    stream.writeStartElement("numberObservations");
    stream.writeCharacters(toString(m_numberObservations));
    stream.writeEndElement();
    stream.writeStartElement("numberImageParameters");
    stream.writeCharacters(toString(m_numberImageParameters));
    stream.writeEndElement();
    stream.writeStartElement("numberConstrainedPointParameters");
    stream.writeCharacters(toString(m_numberConstrainedPointParameters));
    stream.writeEndElement();
    stream.writeStartElement("numberConstrainedImageParameters");
    stream.writeCharacters(toString(m_numberConstrainedImageParameters));
    stream.writeEndElement();
    stream.writeStartElement("numberUnknownParameters");
    stream.writeCharacters(toString(m_numberUnknownParameters));
    stream.writeEndElement();
    stream.writeStartElement("degreesOfFreedom");
    stream.writeCharacters(toString(m_degreesOfFreedom));
    stream.writeEndElement();
    stream.writeStartElement("sigma0");
    stream.writeCharacters(toString(m_sigma0));
    stream.writeEndElement();
    stream.writeStartElement("converged");
    stream.writeCharacters(toString(m_converged));
    stream.writeEndElement();
    stream.writeEndElement();

    stream.writeStartElement("rms");
    stream.writeStartElement("residuals");
    stream.writeAttribute("x", toString(m_rms_rx)); 
    stream.writeAttribute("y", toString(m_rms_ry)); 
    stream.writeAttribute("xy", toString(m_rms_rxy)); 
    stream.writeEndElement(); // end residuals element
    stream.writeStartElement("sigmas");
    stream.writeAttribute("lat", toString(m_rmsSigmaLat)); 
    stream.writeAttribute("lon", toString(m_rmsSigmaLon)); 
    stream.writeAttribute("rad", toString(m_rmsSigmaRad)); 
    stream.writeEndElement(); // end sigmas element
    stream.writeStartElement("imageResidualsLists");
    stream.writeStartElement("residualsList");
    stream.writeAttribute("listSize", toString(m_rmsImageResiduals.size())); 
    for (int i = 0; i < m_rmsImageResiduals.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageResiduals[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageResiduals[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageResiduals[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageResiduals[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageResiduals[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageResiduals[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageResiduals[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageResiduals[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageResiduals[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageResiduals[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageResiduals[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageResiduals[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageResiduals[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageResiduals[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end residuals list
    stream.writeStartElement("sampleResidualsList");
    stream.writeAttribute("listSize", toString(m_rmsImageSampleResiduals.size())); 
    for (int i = 0; i < m_rmsImageSampleResiduals.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageSampleResiduals[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageSampleResiduals[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageSampleResiduals[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageSampleResiduals[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageSampleResiduals[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageSampleResiduals[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageSampleResiduals[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageSampleResiduals[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageSampleResiduals[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageSampleResiduals[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageSampleResiduals[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageSampleResiduals[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageSampleResiduals[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageSampleResiduals[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end sample residuals list
    stream.writeStartElement("lineResidualsList");
    stream.writeAttribute("listSize", toString(m_rmsImageLineResiduals.size())); 
    for (int i = 0; i < m_rmsImageLineResiduals.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageLineResiduals[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageLineResiduals[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageLineResiduals[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageLineResiduals[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageLineResiduals[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageLineResiduals[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageLineResiduals[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageLineResiduals[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageLineResiduals[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageLineResiduals[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageLineResiduals[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageLineResiduals[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageLineResiduals[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageLineResiduals[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end line residuals list
    stream.writeEndElement(); // end image residuals lists
    stream.writeStartElement("imageSigmasLists");
    stream.writeStartElement("xSigmas");
    stream.writeAttribute("listSize", toString(m_rmsImageXSigmas.size())); 
    for (int i = 0; i < m_rmsImageXSigmas.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageXSigmas[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageXSigmas[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageXSigmas[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageXSigmas[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageXSigmas[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageXSigmas[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageXSigmas[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageXSigmas[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageXSigmas[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageXSigmas[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageXSigmas[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageXSigmas[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageXSigmas[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageXSigmas[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end x sigma list
    stream.writeStartElement("ySigmas");
    stream.writeAttribute("listSize", toString(m_rmsImageYSigmas.size())); 
    for (int i = 0; i < m_rmsImageYSigmas.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageYSigmas[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageYSigmas[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageYSigmas[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageYSigmas[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageYSigmas[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageYSigmas[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageYSigmas[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageYSigmas[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageYSigmas[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageYSigmas[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageYSigmas[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageYSigmas[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageYSigmas[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageYSigmas[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end y sigma list
    stream.writeStartElement("zSigmas");
    stream.writeAttribute("listSize", toString(m_rmsImageZSigmas.size())); 
    for (int i = 0; i < m_rmsImageZSigmas.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageZSigmas[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageZSigmas[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageZSigmas[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageZSigmas[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageZSigmas[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageZSigmas[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageZSigmas[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageZSigmas[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageZSigmas[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageZSigmas[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageZSigmas[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageZSigmas[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageZSigmas[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageZSigmas[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end z sigma list
    stream.writeStartElement("raSigmas");
    stream.writeAttribute("listSize", toString(m_rmsImageRASigmas.size())); 
    for (int i = 0; i < m_rmsImageRASigmas.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageRASigmas[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageRASigmas[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageRASigmas[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageRASigmas[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageRASigmas[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageRASigmas[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageRASigmas[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageRASigmas[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageRASigmas[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageRASigmas[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageRASigmas[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageRASigmas[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageRASigmas[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageRASigmas[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end ra sigma list
    stream.writeStartElement("decSigmas");
    stream.writeAttribute("listSize", toString(m_rmsImageDECSigmas.size())); 
    for (int i = 0; i < m_rmsImageDECSigmas.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageDECSigmas[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageDECSigmas[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageDECSigmas[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageDECSigmas[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageDECSigmas[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageDECSigmas[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageDECSigmas[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageDECSigmas[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageDECSigmas[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageDECSigmas[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageDECSigmas[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageDECSigmas[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageDECSigmas[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageDECSigmas[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end dec sigma list
    stream.writeStartElement("twistSigmas");
    stream.writeAttribute("listSize", toString(m_rmsImageTWISTSigmas.size())); 
    for (int i = 0; i < m_rmsImageTWISTSigmas.size(); i++) {
      stream.writeStartElement("statistics");
      stream.writeAttribute("validPixels", toString(m_rmsImageTWISTSigmas[i].ValidPixels()));
      stream.writeAttribute("nullPixels", toString(m_rmsImageTWISTSigmas[i].NullPixels()));
      stream.writeAttribute("lrsPixels", toString(m_rmsImageTWISTSigmas[i].LrsPixels()));
      stream.writeAttribute("lisPixels", toString(m_rmsImageTWISTSigmas[i].LisPixels()));
      stream.writeAttribute("hrsPixels", toString(m_rmsImageTWISTSigmas[i].HrsPixels()));
      stream.writeAttribute("hisPixels", toString(m_rmsImageTWISTSigmas[i].HisPixels()));
      stream.writeAttribute("underRangePixels", toString(m_rmsImageTWISTSigmas[i].UnderRangePixels()));
      stream.writeAttribute("overRangePixels", toString(m_rmsImageTWISTSigmas[i].OverRangePixels()));
      stream.writeAttribute("sum", toString(m_rmsImageTWISTSigmas[i].Sum()));
      stream.writeAttribute("minimum", toString(m_rmsImageTWISTSigmas[i].Minimum()));
      stream.writeAttribute("maximum", toString(m_rmsImageTWISTSigmas[i].Maximum()));
      stream.writeAttribute("validMinimum", toString(m_rmsImageTWISTSigmas[i].ValidMinimum()));
      stream.writeAttribute("validMaximum", toString(m_rmsImageTWISTSigmas[i].ValidMaximum()));
      stream.writeAttribute("removedData", toString(m_rmsImageTWISTSigmas[i].RemovedData()));
      stream.writeEndElement();
    }
    stream.writeEndElement(); // end twist sigma list
    stream.writeEndElement(); // end sigmas lists
    stream.writeEndElement(); // end rms

    stream.writeStartElement("elapsedTime");
    stream.writeAttribute("time", toString(m_elapsedTime)); 
    stream.writeAttribute("errorProp", toString(m_elapsedTimeErrorProp)); 
    stream.writeEndElement();

    stream.writeStartElement("minMaxSigmas");
    stream.writeStartElement("minLat");
    stream.writeAttribute("value", toString(m_minSigmaLatitude)); 
    stream.writeAttribute("pointId", m_minSigmaLatitudePointId); 
    stream.writeEndElement();
    stream.writeStartElement("maxLat");
    stream.writeAttribute("value", toString(m_maxSigmaLatitude)); 
    stream.writeAttribute("pointId", m_maxSigmaLatitudePointId); 
    stream.writeEndElement();
    stream.writeStartElement("minLon");
    stream.writeAttribute("value", toString(m_minSigmaLongitude)); 
    stream.writeAttribute("pointId", m_minSigmaLongitudePointId); 
    stream.writeEndElement();
    stream.writeStartElement("maxLon");
    stream.writeAttribute("value", toString(m_maxSigmaLongitude)); 
    stream.writeAttribute("pointId", m_maxSigmaLongitudePointId); 
    stream.writeEndElement();
    stream.writeStartElement("minRad");
    stream.writeAttribute("value", toString(m_minSigmaRadius)); 
    stream.writeAttribute("pointId", m_minSigmaRadiusPointId); 
    stream.writeEndElement();
    stream.writeStartElement("maxRad");
    stream.writeAttribute("value", toString(m_maxSigmaRadius)); 
    stream.writeAttribute("pointId", m_maxSigmaRadiusPointId); 
    stream.writeEndElement();
    stream.writeEndElement();
    // call max likelihood setup from startElement to fill the rest of these values... 
    stream.writeStartElement("maximumLikelihoodEstimation");
    stream.writeAttribute("numberModels", toString(m_numberMaximumLikelihoodModels)); 
    stream.writeAttribute("maximumLikelihoodIndex", toString(m_maximumLikelihoodIndex)); 
    stream.writeAttribute("maximumLikelihoodMedianR2Residuals", toString(m_maximumLikelihoodMedianR2Residuals)); 
    stream.writeStartElement("cumulativeProbabilityCalculator");
    // TODO ???  m_cumPro->save(stream, project, newProjectRoot);
    stream.writeEndElement();
    stream.writeStartElement("residualsCumulativeProbabilityCalculator");
    // TODO ???  m_cumProRes->save(stream, project, newProjectRoot);
    stream.writeEndElement();
    for (int i = 0; i < m_numberMaximumLikelihoodModels; i++) {
      stream.writeStartElement("model");
      stream.writeAttribute("index", toString(i)); 
      stream.writeAttribute("modelSelection", 
                          MaximumLikelihoodWFunctions::modelToString(m_wFunc[i]->model()));
      stream.writeAttribute("tweakingconstant", toString(m_wFunc[i]->tweakingConstant())); 
      stream.writeAttribute("quantile", toString(m_maximumLikelihoodQuan[i]));
      stream.writeEndElement();
    }
    stream.writeEndElement();
    

  }



  BundleStatistics::XmlHandler::XmlHandler(BundleStatistics *statistics, 
                                                         Project *project) {
    m_statistics = statistics;
    m_project = project;
    m_characters = "";
  }



  BundleStatistics::XmlHandler::~XmlHandler() {
    // ??? compile error ??? delete m_project;
    m_project = NULL;
  }
  


  bool BundleStatistics::XmlHandler::startElement(const QString &namespaceURI, 
                                                                const QString &localName,
                                                                const QString &qName,
                                                                const QXmlAttributes &atts) {
    if (XmlStackedHandler::startElement(namespaceURI, localName, qName, atts)) {
        
      if (localName == "correlationMatrix") {
        // TODO ???  m_statistics->m_correlationMatrix = new CorrelationMatrix(m_project, reader());
      }
      else if (localName == "residuals") {
        
        QString rx = atts.value("x");
        if (!rx.isEmpty()) {
          m_statistics->m_rms_rx = toDouble(rx);
        }

        QString ry = atts.value("y");
        if (!ry.isEmpty()) {
          m_statistics->m_rms_ry = toDouble(ry);
        }

        QString rxy = atts.value("xy");
        if (!rxy.isEmpty()) {
          m_statistics->m_rms_rxy = toDouble(rxy);
        }

      }
      else if (localName == "sigmas") {

        QString lat = atts.value("lat");
        if (!lat.isEmpty()) {
          m_statistics->m_rmsSigmaLat = toDouble(lat);
        }

        QString lon = atts.value("lon");
        if (!lon.isEmpty()) {
          m_statistics->m_rmsSigmaLon = toDouble(lon);
        }

        QString rad = atts.value("rad");
        if (!rad.isEmpty()) {
          m_statistics->m_rmsSigmaRad = toDouble(rad);
        }

      }
      else if (localName == "statisticsResiduals") {

        QString validPixelsStr = atts.value("validPixels");
        QString nullPixelsStr = atts.value("nullPixels");
        QString lrsPixelsStr = atts.value("lrsPixels");
        QString lisPixelsStr = atts.value("lisPixels");
        QString hrsPixelsStr = atts.value("hrsPixels");
        QString hisPixelsStr = atts.value("hisPixels");
        QString underRangePixelsStr = atts.value("underRangePixels");
        QString overRangePixelsStr = atts.value("overRangePixels");
        QString sumStr = atts.value("sum");
        QString minimumStr = atts.value("minimum");
        QString maximumStr = atts.value("maximum");
        QString validMinimumStr = atts.value("validMinimum");
        QString validMaximumStr = atts.value("validMaximum");
        QString removedDataStr = atts.value("removedData");
        if (!validPixelsStr.isEmpty() 
            && !nullPixelsStr.isEmpty()
            && !lrsPixelsStr.isEmpty()
            && !lisPixelsStr.isEmpty()
            && !hrsPixelsStr.isEmpty()
            && !hisPixelsStr.isEmpty()
            && !underRangePixelsStr.isEmpty()
            && !overRangePixelsStr.isEmpty()
            && !sumStr.isEmpty()
            && !minimumStr.isEmpty()
            && !maximumStr.isEmpty()
            && !validMinimumStr.isEmpty()
            && !validMaximumStr.isEmpty()
            && !removedDataStr.isEmpty()
            ) {
          int validPixels = toInt(validPixelsStr);
          int nullPixels = toInt(nullPixelsStr);
          int lrsPixels = toInt(lrsPixelsStr);
          int lisPixels = toInt(lisPixelsStr);
          int hrsPixels = toInt(hrsPixelsStr);
          int hisPixels = toInt(hisPixelsStr);
          int underRangePixels = toInt(underRangePixelsStr);
          int overRangePixels = toInt(overRangePixelsStr);
          double sum = toDouble(sumStr);
          double minimum = toDouble(minimumStr);
          double maximum = toDouble(maximumStr);
          double validMinimum = toDouble(validMinimumStr);
          double validMaximum = toDouble(validMaximumStr);
          bool removedData = toBool(removedDataStr);
          Statistics statsItem;
          statsItem.set(validPixels, nullPixels, lrsPixels, lisPixels, hrsPixels, hisPixels,
                        underRangePixels, overRangePixels, sum, minimum, maximum,
                        validMinimum, validMaximum, removedData);
          m_statistics->m_rmsImageResiduals.append(statsItem);
        }

      }
      else if (localName == "elapsedTime") {

        QString time = atts.value("time");
        if (!time.isEmpty()) {
          m_statistics->m_elapsedTime = toDouble(time);
        }

        QString errorProp = atts.value("errorProp");
        if (!errorProp.isEmpty()) {
          m_statistics->m_elapsedTimeErrorProp = toDouble(errorProp);
        }

      }
      else if (localName == "minLat") {

        QString minLat = atts.value("value");
        if (!minLat.isEmpty()) {
          m_statistics->m_minSigmaLatitude = toDouble(minLat);
        }

        QString minLatPointId = atts.value("pointId");
        if (!minLatPointId.isEmpty()) {
          m_statistics->m_minSigmaLatitudePointId = minLatPointId;
        }

      }
      else if (localName == "maxLat") {

        QString maxLat = atts.value("value");
        if (!maxLat.isEmpty()) {
          m_statistics->m_maxSigmaLatitude = toDouble(maxLat);
        }

        QString maxLatPointId = atts.value("pointId");
        if (!maxLatPointId.isEmpty()) {
          m_statistics->m_maxSigmaLatitudePointId = maxLatPointId;
        }

      }
      else if (localName == "minLon") {

        QString minLon = atts.value("value");
        if (!minLon.isEmpty()) {
          m_statistics->m_minSigmaLongitude = toDouble(minLon);
        }

        QString minLonPointId = atts.value("pointId");
        if (!minLonPointId.isEmpty()) {
          m_statistics->m_minSigmaLongitudePointId = minLonPointId;
        }

      }
      else if (localName == "maxLon") {

        QString maxLon = atts.value("value");
        if (!maxLon.isEmpty()) {
          m_statistics->m_maxSigmaLongitude = toDouble(maxLon);
        }

        QString maxLonPointId = atts.value("pointId");
        if (!maxLonPointId.isEmpty()) {
          m_statistics->m_maxSigmaLongitudePointId = maxLonPointId;
        }

      }
      else if (localName == "minRad") {

        QString minRad = atts.value("value");
        if (!minRad.isEmpty()) {
          m_statistics->m_minSigmaRadius = toDouble(minRad);
        }

        QString minRadPointId = atts.value("pointId");
        if (!minRadPointId.isEmpty()) {
          m_statistics->m_minSigmaRadiusPointId = minRadPointId;
        }

      }
      else if (localName == "maxRad") {

        QString maxRad = atts.value("value");
        if (!maxRad.isEmpty()) {
          m_statistics->m_maxSigmaRadius = toDouble(maxRad);
        }

        QString maxRadPointId = atts.value("pointId");
        if (!maxRadPointId.isEmpty()) {
          m_statistics->m_maxSigmaRadiusPointId = maxRadPointId;
        }

      }
      else if (localName == "maximumLikelihoodEstimation") {

        QString numberModels = atts.value("numberModels");
        if (!numberModels.isEmpty()) {
          m_statistics->m_numberMaximumLikelihoodModels = toInt(numberModels);
        }

        QString maximumLikelihoodIndex = atts.value("maximumLikelihoodIndex");
        if (!maximumLikelihoodIndex.isEmpty()) {
          m_statistics->m_maximumLikelihoodIndex = toInt(maximumLikelihoodIndex);
        }

        QString maximumLikelihoodMedianR2Residuals = atts.value("maximumLikelihoodMedianR2Residuals");
        if (!maximumLikelihoodMedianR2Residuals.isEmpty()) {
          m_statistics->m_maximumLikelihoodMedianR2Residuals = toDouble(maximumLikelihoodMedianR2Residuals);
        }

      }
      else if (localName == "cumulativeProbabilityCalculator") {
        // TODO ???  m_statistics->m_cumPro = new StatCumProbDistDynCalc(m_project, reader());
      }
      else if (localName == "residualsCumulativeProbabilityCalculator") {
        // TODO ???  m_statistics->m_cumProRes = new StatCumProbDistDynCalc(m_project, reader());
      }
      else if (localName == "model") {

        QString index = atts.value("index");
        int i = toInt(index);

        QString model = atts.value("modelSelection");
        QString tweakingConstant = atts.value("tweakingConstant");
        if (!model.isEmpty() && !tweakingConstant.isEmpty()) {
          delete m_statistics->m_wFunc[i];
          m_statistics->m_wFunc[i] = NULL;
          m_statistics->m_wFunc[i] = new MaximumLikelihoodWFunctions(
              MaximumLikelihoodWFunctions::stringToModel(model),
              toDouble(tweakingConstant));
        }

        QString quantile = atts.value("quantile");
        if (!quantile.isEmpty()) {
          m_statistics->m_maximumLikelihoodQuan[i] = toDouble(quantile);
          // m_bundleSettings->m_maximumLikelihood.append(
          //    qMakePair(MaximumLikelihoodWFunctions::stringToModel(type),
          //              toDouble(quantile)));
        }

      }
    }
    return true;
  }



  bool BundleStatistics::XmlHandler::characters(const QString &ch) {
    m_characters += ch;
    return XmlStackedHandler::characters(ch);
  }



  bool BundleStatistics::XmlHandler::endElement(const QString &namespaceURI, const QString &localName,
                                     const QString &qName) {
    if (!m_characters.isEmpty()) {
      if (localName == "id") {
        delete m_statistics->m_id;
        m_statistics->m_id = NULL;
        m_statistics->m_id = new QUuid(m_characters);
      }
      if (localName == "instrumentId") {
        m_statistics->m_instrumentId = m_characters;
      }
      if (localName == "numberFixedPoints") {
        m_statistics->m_numberFixedPoints = toInt(m_characters);
      }
      if (localName == "numberIgnoredPoints") {
        m_statistics->m_numberIgnoredPoints = toInt(m_characters);
      }
      if (localName == "numberHeldImages") {
        m_statistics->m_numberHeldImages = toInt(m_characters);
      }
      if (localName == "rejectionLimit") {
        m_statistics->m_rejectionLimit = toDouble(m_characters);
      }
      if (localName == "numberRejectedObservations") {
        m_statistics->m_numberRejectedObservations = toInt(m_characters);
      }
      if (localName == "numberObservations") {
        m_statistics->m_numberObservations = toInt(m_characters);
      }
      if (localName == "numberImageParameters") {
        m_statistics->m_numberImageParameters = toInt(m_characters);
      }
      if (localName == "numberConstrainedPointParameters") {
        m_statistics->m_numberConstrainedPointParameters = toInt(m_characters);
      }
      if (localName == "numberConstrainedImageParameters") {
        m_statistics->m_numberConstrainedImageParameters = toInt(m_characters);
      }
      if (localName == "numberUnknownParameters") {
        m_statistics->m_numberUnknownParameters = toInt(m_characters);
      }
      if (localName == "degreesOfFreedom") {
        m_statistics->m_degreesOfFreedom = toInt(m_characters);
      }
      if (localName == "sigma0") {
        m_statistics->m_sigma0 = toDouble(m_characters);
      }
      if (localName == "converged") {
        m_statistics->m_converged = toBool(m_characters);
      }
      m_characters = "";
    }
    return XmlStackedHandler::endElement(namespaceURI, localName, qName);
  }



  QDataStream &BundleStatistics::write(QDataStream &stream) const {
    stream << m_id->toString()
           << *m_correlationMatrix
           << (qint32)m_numberFixedPoints
           << (qint32)m_numberIgnoredPoints
           << (qint32)m_numberHeldImages
           << m_rms_rx << m_rms_ry << m_rms_rxy
           << m_rejectionLimit
           << (qint32)m_numberRejectedObservations
           << (qint32)m_numberObservations
           << (qint32)m_numberImageParameters
           << (qint32)m_numberConstrainedPointParameters
           << (qint32)m_numberConstrainedImageParameters
           << (qint32)m_numberUnknownParameters
           << (qint32)m_degreesOfFreedom
           << m_sigma0
           << m_elapsedTime << m_elapsedTimeErrorProp
           << m_converged
           << m_rmsImageSampleResiduals << m_rmsImageLineResiduals
           << m_rmsImageResiduals
           << m_rmsImageXSigmas << m_rmsImageYSigmas << m_rmsImageZSigmas
           << m_rmsImageRASigmas << m_rmsImageDECSigmas << m_rmsImageTWISTSigmas
           << m_minSigmaLatitude << m_minSigmaLatitudePointId
           << m_maxSigmaLatitude << m_maxSigmaLatitudePointId
           << m_minSigmaLongitude << m_minSigmaLongitudePointId
           << m_maxSigmaLongitude << m_maxSigmaLongitudePointId
           << m_minSigmaRadius << m_minSigmaRadiusPointId
           << m_maxSigmaRadius << m_maxSigmaRadiusPointId
           << m_rmsSigmaLat << m_rmsSigmaLon << m_rmsSigmaRad
           << (qint32)m_numberMaximumLikelihoodModels << *m_wFunc[0] << *m_wFunc[1] << *m_wFunc[2]
           << m_maximumLikelihoodQuan[0] << m_maximumLikelihoodQuan[1] << m_maximumLikelihoodQuan[2]
           << (qint32)m_maximumLikelihoodIndex << *m_cumPro << *m_cumProRes
           << m_maximumLikelihoodMedianR2Residuals;
    return stream;
  }



  QDataStream &BundleStatistics::read(QDataStream &stream) {
    QString id;
    qint32 numberFixedPoints, numberIgnoredPoints, numberHeldImages, numberRejectedObservations,
           numberObservations, numberImageParameters, numberConstrainedPointParameters,
           numberConstrainedImageParameters, numberUnknownParameters, degreesOfFreedom,
           numberMaximumLikelihoodModels, maximumLikelihoodIndex;
    CorrelationMatrix correlationMatrix;
    QVector<MaximumLikelihoodWFunctions> wFunc;
    wFunc.resize(3);
    StatCumProbDistDynCalc cumPro;
    StatCumProbDistDynCalc cumProRes;
    stream >> id 
           >> correlationMatrix
           >> numberFixedPoints
           >> numberIgnoredPoints
           >> numberHeldImages
           >> m_rms_rx >> m_rms_ry >> m_rms_rxy
           >> m_rejectionLimit
           >> numberRejectedObservations
           >> numberObservations
           >> numberImageParameters
           >> numberConstrainedPointParameters
           >> numberConstrainedImageParameters
           >> numberUnknownParameters
           >> degreesOfFreedom
           >> m_sigma0
           >> m_elapsedTime >> m_elapsedTimeErrorProp
           >> m_converged
           >> m_rmsImageSampleResiduals >> m_rmsImageLineResiduals
           >> m_rmsImageResiduals
           >> m_rmsImageXSigmas >> m_rmsImageYSigmas >> m_rmsImageZSigmas
           >> m_rmsImageRASigmas >> m_rmsImageDECSigmas >> m_rmsImageTWISTSigmas
           >> m_minSigmaLatitude >> m_minSigmaLatitudePointId
           >> m_maxSigmaLatitude >> m_maxSigmaLatitudePointId
           >> m_minSigmaLongitude >> m_minSigmaLongitudePointId
           >> m_maxSigmaLongitude >> m_maxSigmaLongitudePointId
           >> m_minSigmaRadius >> m_minSigmaRadiusPointId
           >> m_maxSigmaRadius >> m_maxSigmaRadiusPointId
           >> m_rmsSigmaLat >> m_rmsSigmaLon >> m_rmsSigmaRad
           >> numberMaximumLikelihoodModels >> wFunc[0] >> wFunc[1] >> wFunc[2]
           >> m_maximumLikelihoodQuan[0] >> m_maximumLikelihoodQuan[1] >> m_maximumLikelihoodQuan[2]
           >> maximumLikelihoodIndex
           >> cumPro >> cumProRes
           >> m_maximumLikelihoodMedianR2Residuals;

    m_id = new QUuid(id);
    m_numberFixedPoints                = (int)numberFixedPoints;
    m_numberIgnoredPoints              = (int)numberIgnoredPoints;
    m_numberHeldImages                 = (int)numberHeldImages;
    m_numberRejectedObservations       = (int)numberRejectedObservations;
    m_numberObservations               = (int)numberObservations;
    m_numberImageParameters            = (int)numberImageParameters;
    m_numberConstrainedPointParameters = (int)numberConstrainedPointParameters;
    m_numberConstrainedImageParameters = (int)numberConstrainedImageParameters;
    m_numberUnknownParameters          = (int)numberUnknownParameters;
    m_degreesOfFreedom                 = (int)degreesOfFreedom;
    m_numberMaximumLikelihoodModels    = (int)numberMaximumLikelihoodModels;
    m_maximumLikelihoodIndex           = (int)maximumLikelihoodIndex;

    delete m_correlationMatrix;
    m_correlationMatrix = new CorrelationMatrix(correlationMatrix);

    for (int i = 0; i < m_numberMaximumLikelihoodModels; i++){

      delete m_wFunc[i];
      m_wFunc[i] = new MaximumLikelihoodWFunctions(wFunc[i]);
    }

    delete m_cumPro;
    m_cumPro = new StatCumProbDistDynCalc(cumPro);

    delete m_cumProRes;
    m_cumProRes = new StatCumProbDistDynCalc(cumProRes);
    
    return stream;
  }



  QDataStream &operator<<(QDataStream &stream, const BundleStatistics &bundleStatistics) {
    return bundleStatistics.write(stream);
  }



  QDataStream &operator>>(QDataStream &stream, BundleStatistics &bundleStatistics) {
    return bundleStatistics.read(stream);
  }   



  /**
   * Writes matrix to binary disk file pointed to by QDataStream stream
   */



  // ??? QList< QPair< MaximumLikelihoodWFunctions, double > >
  // ???     BundleStatistics::maximumLikelihoodModels() const {
  // ??? }

}





#if 0
void BundleStatistics::setNumberHeldImages(SerialNumberList m_heldSnList,
                                           SerialNumberList *snList) {
  for (int i = 0; i < snList->Size(); i++) {
    if (m_heldSnList->HasSerialNumber(snList->SerialNumber(i))) {
      m_heldImages++;
    }
  }

}






/**
 * This method computes the focal plane residuals for the measures.
 *
 * @history 2012-01-18 Debbie A. Cook - Fixed the computation of vx
 *                            and vy to make sure they are focal
 *                            plane x and y residuals instead of 
 *                            image sample and line residuals.
 *
 */
double BundleStatistics::computeResiduals(
                                          ControlNet *cnet,
                                          std::vector< boost::numeric::ublas::bounded_vector< double, 3 > > pointWeights,
                                          std::vector< boost::numeric::ublas::bounded_vector< double, 3 > > pointCorrections,
                                          boost::numeric::ublas::vector< double >imageCorrections,
                                          std::vector< double > imageParameterWeights,
                                          int numImagePartials,
                                          int rank) {
  double vtpv         = 0.0;
  double vtpv_control = 0.0;
  double vtpv_image   = 0.0;
  double dWeight;
  double v, vx, vy;

  // clear residual stats vectors
  m_statsrx.Reset();
  m_statsry.Reset();
  m_statsrxy.Reset();

  // vtpv for image coordinates
  int numObjectPoints = cnet->GetNumPoints();

  for (int i = 0; i < numObjectPoints; i++) {
    ControlPoint *point = cnet->GetPoint(i);
    if (point->IsIgnored()) {
      continue;
    }

    point->ComputeResiduals();

    int numMeasures = point->GetNumMeasures();
    for (int j = 0; j < numMeasures; j++) {
      const ControlMeasure *measure = point->GetMeasure(j);
      if (measure->IsIgnored()) {
        continue;
      }

      dWeight = 1.4 * (measure->Camera())->PixelPitch();
      dWeight = 1.0 / dWeight;
      dWeight *= dWeight;

      vx = measure->GetFocalPlaneMeasuredX() - measure->GetFocalPlaneComputedX();
      vy = measure->GetFocalPlaneMeasuredY() - measure->GetFocalPlaneComputedY();

      // if rejected, don't include in statistics
      if (measure->IsRejected()) {
        continue;
      }

      m_statsrx.AddData(vx);
      m_statsry.AddData(vy);
      m_statsrxy.AddData(vx);
      m_statsrxy.AddData(vy);


      vtpv += vx * vx * dWeight + vy * vy * dWeight;
    }
  }


  // add vtpv from constrained 3D points
  int nPointIndex = 0;
  for (int i = 0; i < numObjectPoints; i++) {
    const ControlPoint *point = cnet->GetPoint(i);
    if (point->IsIgnored()) {
      continue;
    }

    // get weight and correction vector for this point
    bounded_vector< double, 3 > &weights = pointWeights[nPointIndex];
    bounded_vector< double, 3 > &corrections = pointCorrections[nPointIndex];

    if (weights[0] > 0.0) {
      vtpv_control += corrections[0] * corrections[0] * weights[0];
    }
    if (weights[1] > 0.0) {
      vtpv_control += corrections[1] * corrections[1] * weights[1];
    }
    if (weights[2] > 0.0) {
      vtpv_control += corrections[2] * corrections[2] * weights[2];
    }

    nPointIndex++;
  }

  // add vtpv from constrained image parameters
  int n = 0;
  do {
    for (int j = 0; j < numImagePartials; j++) {
      if (imageParameterWeights[j] > 0.0) {
        v = imageCorrections[n];
        vtpv_image += v * v * imageParameterWeights[j];
      }

      n++;
    }

  }
  while (n < rank);

  vtpv = vtpv + vtpv_control + vtpv_image;

  // Compute rms for all image coordinate residuals
  // separately for x, y, then x and y together
  m_rms_rx =  m_statsrx.Rms();
  m_rms_ry =  m_statsry.Rms();
  m_rms_rxy = m_statsrxy.Rms();

  return vtpv;
}

/**
* Compute rejection limit.
*
*/
double BundleStatistics::computeRejectionLimit(ControlNet *cnet,
                                               double outlierRejectionMultiplier,
                                               int numObservations) {

  double vx, vy;
  m_numberObservations = numObservations;
  int nResiduals = m_numberObservations / 2; // make sure this has been set???

  QVector< double > resvectors;
  resvectors.resize(nResiduals);

  // load magnitude (squared) of residual vector
  int numObs = 0;
  int numObjectPoints = cnet->GetNumPoints();
  for (int i = 0; i < numObjectPoints; i++) {

    const ControlPoint *point = cnet->GetPoint(i);
    if (point->IsIgnored()) {
      continue;
    }

    if (point->IsRejected()) {
      continue;
    }

    int numMeasures = point->GetNumMeasures();
    for (int j = 0; j < numMeasures; j++) {

      const ControlMeasure *measure = point->GetMeasure(j);
      if (measure->IsIgnored()) {
        continue;
      }

      if (measure->IsRejected()) {
        continue;
      }

      vx = measure->GetSampleResidual();
      vy = measure->GetLineResidual();

      resvectors[numObs] = sqrt(vx * vx + vy * vy);

      numObs++;
    }
  }

  // sort vectors
  std::sort(resvectors.begin(), resvectors.end());

  double median;
  double mediandev;
  double mad;

  int    nmidpoint = nResiduals / 2;

  if (nResiduals % 2 == 0) {
    median = (resvectors[nmidpoint - 1] + resvectors[nmidpoint]) / 2;
  }
  else {
    median = resvectors[nmidpoint];
  }

  // compute M.A.D.
  for (int i = 0; i < nResiduals; i++) {
    resvectors[i] = fabs(resvectors[i] - median);
  }

  std::sort(resvectors.begin(), resvectors.end());

  if (nResiduals % 2 == 0) {
    mediandev = (resvectors[nmidpoint - 1] + resvectors[nmidpoint]) / 2;
  }
  else {
    mediandev = resvectors[nmidpoint];
  }

  mad = 1.4826 * mediandev;

  return (median + outlierRejectionMultiplier * mad);
}



/**
 * Flag outlier measurements.
 *
 */
bool BundleStatistics::flagOutliers(ControlNet *cnet) {
  double vx, vy;
  int    nRejected;
  int    ntotalrejected      = 0;

  int    nIndexMaxResidual;
  double dMaxResidual;
  double dSumSquares;
  double dUsedRejectionLimit = m_rejectionLimit;

//    if ( m_rejectionLimit < 0.05 )
//        dUsedRejectionLimit = 0.14;
//        dUsedRejectionLimit = 0.05;

  int    nComingBack         = 0;

  int    numObjectPoints       = cnet->GetNumPoints();
  for (int i = 0; i < numObjectPoints; i++) {
    ControlPoint *point = cnet->GetPoint(i);
    if (point->IsIgnored()) {
      continue;
    }

    point->ZeroNumberOfRejectedMeasures();

    nRejected = 0;
    nIndexMaxResidual = -1;
    dMaxResidual = -1.0;

    int numMeasures = point->GetNumMeasures();
    for (int j = 0; j < numMeasures; j++) {

      ControlMeasure *measure = point->GetMeasure(j);
      if (measure->IsIgnored()) {
        continue;
      }

      vx = measure->GetSampleResidual();
      vy = measure->GetLineResidual();

      dSumSquares = sqrt(vx * vx + vy * vy);

      // measure is good
      if (dSumSquares <= dUsedRejectionLimit) {

        // was it previously rejected?
        if (measure->IsRejected()) {
          printf("Coming back in: %s\r", point->GetId().toAscii().data());
          nComingBack++;
          cnet->DecrementNumberOfRejectedMeasuresInImage(measure->GetCubeSerialNumber());
        }

        measure->SetRejected(false);
        continue;
      }

      // if it's still rejected, skip it
      if (measure->IsRejected()) {
        nRejected++;
        ntotalrejected++;
        continue;
      }

      if (dSumSquares > dMaxResidual) {
        dMaxResidual = dSumSquares;
        nIndexMaxResidual = j;
      }
    }

    // no observations above the current rejection limit for this 3D point
    if (dMaxResidual == -1.0 || dMaxResidual <= dUsedRejectionLimit) {
      point->SetNumberOfRejectedMeasures(nRejected);
      continue;
    }

    // this is another kluge - if we only have two observations
    // we won't reject (for now)
    if ((numMeasures - (nRejected + 1)) < 2) {
      point->SetNumberOfRejectedMeasures(nRejected);
      continue;
    }

    // otherwise, we have at least one observation
    // for this point whose residual is above the
    // current rejection limit - we'll flag the
    // worst of these as rejected
    ControlMeasure *rejected = point->GetMeasure(nIndexMaxResidual);
    rejected->SetRejected(true);
    nRejected++;
    point->SetNumberOfRejectedMeasures(nRejected);
    cnet->IncrementNumberOfRejectedMeasuresInImage(rejected->GetCubeSerialNumber());
    ntotalrejected++;

    // do we still have sufficient remaining observations for this 3D point?
    if ((numMeasures - nRejected) < 2) {
      point->SetRejected(true);
      printf("Rejecting Entire Point: %s\r", point->GetId().toAscii().data());
    }
    else point->SetRejected(false);

//      int ndummy = point->GetNumberOfRejectedMeasures();
//      printf("Rejected for point %s = %d\n", point->GetId().toAscii().data(), ndummy);
//      printf("%s: %20.10lf  %20.10lf*\n",point->GetId().toAscii().data(), rejected->GetSampleResidual(), rejected->GetLineResidual());
  }

  m_numberRejectedObservations = 2 * ntotalrejected;

  printf("\n\t       Rejected Observations:%10d (Rejection Limit:%12.5lf\n",
         m_numberRejectedObservations, dUsedRejectionLimit);

//    std::cout << "Measures that came back: " << nComingBack << std::endl;

  return true;
}


#endif

#if 0
void setNumberImageParameters(int numImagePartials, int observations) {
  m_ImageParameters = numImagePartials * observations;
}
#endif



