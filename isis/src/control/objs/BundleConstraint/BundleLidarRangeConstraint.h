#ifndef BundleLidarRangeConstraint_h
#define BundleLidarRangeConstraint_h
/**
 * @file
 * $Revision: 1.17 $
 * $Date: 2010/03/27 07:01:33 $
 *
 *   Unless noted otherwise, the portions of Isis written by the USGS are public
 *   domain. See individual third-party library and package descriptions for
 *   intellectual property information,user agreements, and related information.
 *
 *   Although Isis has been used by the USGS, no warranty, expressed or implied,
 *   is made by the USGS as to the accuracy and functioning of such software
 *   and related material nor shall the fact of distribution constitute any such
 *   warranty, and no responsibility is assumed by the USGS in connection
 *   therewith.
 *
 *   For additional information, launch
 *   $ISISROOT/doc//documents/Disclaimers/Disclaimers.html in a browser or see
 *   the Privacy &amp; Disclaimers page on the Isis website,
 *   http://isis.astrogeology.usgs.gov, and the USGS privacy and disclaimers on
 *   http://www.usgs.gov/privacy.html.
 */
// Qt Library
#include <QSharedPointer>

// Isis Library
#include "BundleConstraint.h"
#include "LinearAlgebra.h"
#include "SparseBlockMatrix.h"

namespace Isis {
  class BundleMeasure;
  typedef QSharedPointer<BundleMeasure> BundleMeasureQsp;
  class BundleObservation;
  typedef QSharedPointer<BundleObservation> BundleObservationQsp;
  class BundleLidarControlPoint;
  typedef QSharedPointer<BundleLidarControlPoint> BundleLidarControlPointQsp;

  /**
   * @brief Implements range constraint between lidar surface point and simultaneous image for
   * bundle adjustment.
   *
   * @ingroup Control
   *
   * @author 2018-03-05 Ken Edmundson
   *
   * @internal
   *   @history 2018-03-05 Ken Edmundson - Original version.
   */
  class BundleLidarRangeConstraint : public BundleConstraint {
    public:
      // constructors
      BundleLidarRangeConstraint();
      BundleLidarRangeConstraint(BundleMeasureQsp bundleMeasure);

      // copy constructor
      BundleLidarRangeConstraint(const BundleLidarRangeConstraint &src);

      // destructor
      ~BundleLidarRangeConstraint();

      // Assignment operator
      BundleLidarRangeConstraint &operator= (const BundleLidarRangeConstraint &src);

      bool formConstraint(LinearAlgebra::MatrixUpperTriangular &N22,
                          SparseBlockColumnMatrix &N12,
                          LinearAlgebra::Vector &n2,
                          LinearAlgebra::VectorCompressed &n1,
                          SparseBlockMatrix& sparseNormals);

      void updateRightHandSide();
      SparseBlockMatrix &normalsSpkMatrix();
      LinearAlgebra::Vector &rightHandSideVector();

      QString formatBundleOutputString();

    private:
      void constructMatrices();
      void positionContinuity(int &designRow);

      //! parent BundleLidarControlPoint
      BundleMeasureQsp m_bundleMeasure;
      BundleObservationQsp m_bundleObservation;
      BundleLidarControlPoint* m_bundleLidarControlPoint;

      int m_numberCoefficients;                                //! # coefficients

      LinearAlgebra::MatrixCompressed m_designMatrix;          //! design matrix
      SparseBlockMatrix m_normalsSpkMatrix;                    //! normals contribution to position
      LinearAlgebra::Vector m_rightHandSide;                   //! right hand side of normals
      LinearAlgebra::Vector m_omcVector;                       //! observed minus corrected vector

      LinearAlgebra::Matrix m_coeffRangeImage;
      LinearAlgebra::Matrix m_coeffRangePoint;
      LinearAlgebra::Vector m_coeffRangeRHS;
  };

  //! Typdef for BundleLidarRangeConstraint QSharedPointer.
  typedef QSharedPointer<BundleLidarRangeConstraint> BundleLidarRangeConstraintQsp;
};

#endif
