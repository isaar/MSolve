using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Materials.Crack
{
    class HomogeneousIsotropicMaterialFactory2D: ICrackMaterialFactory2D
    {
        private readonly double poissonRatio;
        private readonly double shearModulus;
        private readonly double criticalFractureTougnhess;

        public HomogeneousIsotropicMaterialFactory2D(double youngModulus, double poissonRatio, 
            double criticalFractureTougnhess)
        {
            if (youngModulus < 0.0)
            {
                throw new ArgumentException("Young's modulus must be non negative but was: " + youngModulus);
            }
            if (poissonRatio < -1.0 || poissonRatio > 0.5)
            {
                throw new ArgumentException("Poisson's ratio must be in the range [-1, 0.5] but was: " + poissonRatio);
            }
            if (criticalFractureTougnhess <= 0.0)
            {
                throw new ArgumentException("The critical fracture toughness must be positive but was: "
                    + criticalFractureTougnhess);
            }

            this.poissonRatio = poissonRatio;
            this.shearModulus = 0.5 * youngModulus / (1.0 + poissonRatio);
            this.criticalFractureTougnhess = criticalFractureTougnhess;
        }

        public CrackMaterial2D FindMaterialAtPoint(ICartesianPoint2D point)
        {
            return new CrackMaterial2D(poissonRatio, shearModulus, criticalFractureTougnhess);
        }
    }
}
