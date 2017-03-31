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
        public enum PlaneProblem
        {
            PlaneStrain, PlaneStress  
        }

        private readonly double poissonRatio;
        private readonly double youngModulus;
        private readonly double criticalFractureTougnhess;
        private readonly double kolosovCoefficient;

        public HomogeneousIsotropicMaterialFactory2D(double youngModulus, double poissonRatio, 
            double criticalFractureTougnhess, PlaneProblem planeProblem)
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

            this.youngModulus = youngModulus;
            this.poissonRatio = poissonRatio;
            this.criticalFractureTougnhess = criticalFractureTougnhess;

            if (planeProblem == PlaneProblem.PlaneStrain) kolosovCoefficient = 3 - 4 * poissonRatio;
            else kolosovCoefficient = (3 - poissonRatio) / (1 + poissonRatio);
        }

        public CrackMaterial2D FindMaterialAtPoint(ICartesianPoint2D point)
        {
            return new CrackMaterial2D(youngModulus, poissonRatio, 
                0.5 * youngModulus / (1.0 + poissonRatio), criticalFractureTougnhess, kolosovCoefficient);
        }
    }
}
