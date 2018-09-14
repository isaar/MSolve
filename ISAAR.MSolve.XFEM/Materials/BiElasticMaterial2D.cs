using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Enrichments.Items;


namespace ISAAR.MSolve.XFEM.Materials
{
    class BiElasticMaterial2D: IMaterialField2D
    {
        private readonly MaterialInterface2D bimaterialInterface;

        // TODO: Perhaps these can be stored as a MaterialPoint2D and then return it (instead of copying to new)
        public double YoungModulus1 { get; }
        public double EquivalentYoungModulus1 { get; }
        public double PoissonRatio1 { get; }
        public double EquivalentPoissonRatio1 { get; }
        public double Thickness1 { get; }

        public double YoungModulus2 { get; }
        public double EquivalentYoungModulus2 { get; }
        public double PoissonRatio2 { get; }
        public double EquivalentPoissonRatio2 { get; }
        public double Thickness2 { get; }

        /// <summary>
        /// Material properties for the "positive" subdomain designed by the interface are denoted as 1.
        /// Material properties for the "negative" subdomain designed by the interface are denoted as 2. 
        /// </summary>
        /// <param name="youngModulus1"></param>
        /// <param name="poissonRatio1"></param>
        /// <param name="youngModulus2"></param>
        /// <param name="poissonRatio2"></param>
        /// <param name="bimaterialInterface"></param>
        /// <returns></returns>
        public static BiElasticMaterial2D CreateMaterialForPlainStrain(double youngModulus1, double poissonRatio1,
             double youngModulus2, double poissonRatio2, MaterialInterface2D bimaterialInterface)
        {
            double equivalentE1 = youngModulus1 / (1.0 - poissonRatio1 * poissonRatio1);
            double equivalentV1 = poissonRatio1 / (1.0 - poissonRatio1);
            double thickness1 = 1.0;
            double equivalentE2 = youngModulus2 / (1.0 - poissonRatio2 * poissonRatio2);
            double equivalentV2 = poissonRatio2 / (1.0 - poissonRatio2);
            double thickness2 = 1.0;
            return new BiElasticMaterial2D(youngModulus1, equivalentE1, poissonRatio1, equivalentV1, thickness1,
                youngModulus2, equivalentE2, poissonRatio2, equivalentV2, thickness2, bimaterialInterface);
        }

        /// <summary>
        /// Material properties for the "positive" subdomain designed by the interface are denoted as 1.
        /// Material properties for the "negative" subdomain designed by the interface are denoted as 2. 
        /// </summary>
        /// <param name="youngModulus1"></param>
        /// <param name="poissonRatio1"></param>
        /// <param name="youngModulus2"></param>
        /// <param name="poissonRatio2"></param>
        /// <param name="bimaterialInterface"></param>
        public static BiElasticMaterial2D CreateMaterialForPlainStress(double youngModulus1, double poissonRatio1, 
            double thickness1, double youngModulus2, double poissonRatio2, double thickness2, 
            MaterialInterface2D bimaterialInterface)
        {
            return new BiElasticMaterial2D(youngModulus1, youngModulus1, poissonRatio1, poissonRatio1, thickness1,
                youngModulus2, youngModulus2, poissonRatio2, poissonRatio2, thickness2, bimaterialInterface);
        }

        private BiElasticMaterial2D(double youngModulus1, double equivalentYoungModulus1, double poissonRatio1,
            double equivalentPoissonRatio1, double thickness1, double youngModulus2, double equivalentYoungModulus2,
            double poissonRatio2, double equivalentPoissonRatio2, double thickness2,
            MaterialInterface2D bimaterialInterface)
        {
            MaterialUtilities.CheckYoungModulus(youngModulus1);
            MaterialUtilities.CheckPoissonRatio(poissonRatio1);
            MaterialUtilities.CheckThickness(thickness1);
            MaterialUtilities.CheckYoungModulus(youngModulus2);
            MaterialUtilities.CheckPoissonRatio(poissonRatio2);
            MaterialUtilities.CheckThickness(thickness2);

            this.YoungModulus1 = youngModulus1;
            EquivalentYoungModulus1 = equivalentYoungModulus1;
            this.PoissonRatio1 = poissonRatio1;
            this.EquivalentPoissonRatio1 = equivalentPoissonRatio1;
            this.Thickness1 = thickness1;

            this.YoungModulus2 = youngModulus2;
            EquivalentYoungModulus2 = equivalentYoungModulus2;
            this.PoissonRatio2 = poissonRatio2;
            this.EquivalentPoissonRatio2 = equivalentPoissonRatio2;
            this.Thickness2 = thickness2;

            this.bimaterialInterface = bimaterialInterface;
        }

        public double GetYoungModulusAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        {
            return IsMaterial1(point, interpolation) ? YoungModulus1 : YoungModulus2;
        }

        public double GetEquivalentYoungModulusAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        {
            return IsMaterial1(point, interpolation) ? EquivalentYoungModulus1 : EquivalentYoungModulus2;
        }

        public double GetPoissonRatioAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        {
            return IsMaterial1(point, interpolation) ? PoissonRatio1 : PoissonRatio2;
        }

        public double GetEquivalentPoissonRatioAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        {
            return IsMaterial1(point, interpolation) ? EquivalentPoissonRatio1 : EquivalentPoissonRatio2;
        }

        public double GetThicknessAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        {
            return IsMaterial1(point, interpolation) ? Thickness1 : Thickness2;
        }

        public Matrix CalculateConstitutiveMatrixAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        {
            double eqE, eqV;
            if (IsMaterial1(point, interpolation))
            {
                eqE = EquivalentYoungModulus1;
                eqV = EquivalentPoissonRatio1;
            }
            else
            {
                eqE = EquivalentYoungModulus2;
                eqV = EquivalentPoissonRatio2;
            }

            double scalar = eqE / (1 - eqV * eqV);
            var matrix = Matrix.CreateZero(3, 3);
            matrix[0, 0] = scalar;
            matrix[0, 1] = scalar * eqV;
            matrix[1, 0] = scalar * eqV;
            matrix[1, 1] = scalar;
            matrix[2, 2] = 0.5 * eqE / (1 + eqV);
            return matrix;
        }

        public void UpdateDistributions() { } // Do nothing for homogeneous preperties
        public void UpdateStrains() { } // Do nothing for elastic properties.

        private bool IsMaterial1(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        {
            ICartesianPoint2D cartesianPoint = interpolation.TransformPointNaturalToGlobalCartesian(point);
            MaterialInterface2D.Subdomain subdomain = bimaterialInterface.LocatePoint(cartesianPoint);
            if (subdomain == MaterialInterface2D.Subdomain.Positive) return true;
            else if (subdomain == MaterialInterface2D.Subdomain.Negative) return false;
            else throw new ArgumentException("The point (xi, eta) = " + point + " - (x, y) = " + cartesianPoint + 
                ", lies on the bi-material interface");
        }
    }
}
