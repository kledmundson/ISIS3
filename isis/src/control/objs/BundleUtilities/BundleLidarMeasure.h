#ifndef BundleLidarMeasure_h
#define BundleLidarMeasure_h
/**
 * @file
 * $Revision: 1.20 $
 * $Date: 2014/5/22 01:35:17 $
 *
 *   Unless noted otherwise, the portions of Isis written by the USGS are
 *   public domain. See individual third-party library and package descriptions
 *   for intellectual property information, user agreements, and related
 *   information.
 *
 *   Although Isis has been used by the USGS, no warranty, expressed or
 *   implied, is made by the USGS as to the accuracy and functioning of such
 *   software and related material nor shall the fact of distribution
 *   constitute any such warranty, and no responsibility is assumed by the
 *   USGS in connection therewith.
 *
 *   For additional information, launch
 *   $ISISROOT/doc//documents/Disclaimers/Disclaimers.html
 *   in a browser or see the Privacy &amp; Disclaimers page on the Isis website,
 *   http://isis.astrogeology.usgs.gov, and the USGS privacy and disclaimers on
 *   http://www.usgs.gov/privacy.html.
 */
// Qt Library
#include <QSharedPointer>

// Isis Library
#include "BundleMeasure.h"
#include "LinearAlgebra.h"
#include "SparseBlockMatrix.h"


namespace Isis {
  class BundleControlPoint;
  class BundleLidarControlPoint;
  class ControlMeasure;

/**
   * @brief A container class for a LidarControlMeasure.
   *
   * TODO: update description
   *
   * This class is used as a wrapper around a ControlMeasure to provide the necessary information
   * for BundleAdjust. This class can be used to get the parent bundle observation solve settings 
   * for observation mode adjustment.
   *
   * Note that a BundleLidarMeasure should be created from a non-ignored ControlMeasure.
   *
   * @ingroup ControlNetworks
   *
   * @author 2018-03-07 Ken Edmundson
   *
   * @internal
   *   @history 2018-03-07 Ken Edmundson - Original version.
   */

  class BundleLidarMeasure: public BundleMeasure {

    public:
      // constructor
      BundleLidarMeasure(ControlMeasure *controlMeasure,
                         BundleLidarControlPoint *bundleControlPoint);

      // copy constructor
//      BundleLidarMeasure(const BundleLidarMeasure &src);

      // destructor
      ~BundleLidarMeasure();

      BundleLidarMeasure &operator=(const BundleLidarMeasure &src);

      bool applyLidarRangeConstraint(LinearAlgebra::MatrixUpperTriangular &N22,
                                     SparseBlockColumnMatrix &N12,
                                     LinearAlgebra::Vector &n2,
                                     LinearAlgebra::VectorCompressed &n1,
                                     SparseBlockMatrix &sparseNormals);

    private:
  };
  //! Definition for BundleLidarMeasureQsp, a shared pointer to a BundleLidarMeasure.
  typedef QSharedPointer<BundleLidarMeasure> BundleLidarMeasureQsp;
}

#endif // BundleLidarMeasure_h
