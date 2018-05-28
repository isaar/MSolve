using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;

namespace ISAAR.MSolve.XFEM.Materials
{
    class HomogeneousElasticMaterial2D: IMaterialField2D
    {
        public double HomogeneousYoungModulus { get; }
        public double HomogeneousEquivalentYoungModulus { get; }
        public double HomogeneousPoissonRatio { get; }
        public double HomogeneousEquivalentPoissonRatio { get; }
        public double HomogeneousThickness { get; }

        public static HomogeneousElasticMaterial2D CreateMaterialForPlaneStrain(double youngModulus, 
            double poissonRatio)
        {
            double equivalentE = youngModulus / (1.0 - poissonRatio * poissonRatio);
            double equivalentV = poissonRatio / (1.0 - poissonRatio);
            double thickness = 1.0;
            return new HomogeneousElasticMaterial2D(youngModulus, equivalentE, poissonRatio, equivalentV, thickness);
        }

        public static HomogeneousElasticMaterial2D CreateMaterialForPlaneStress(double youngModulus,
            double poissonRatio, double thickness)
        {
            return new HomogeneousElasticMaterial2D(youngModulus, youngModulus, poissonRatio, poissonRatio, thickness);
        }

        private HomogeneousElasticMaterial2D(double youngModulus, double equivalentYoungModulus, double poissonRatio,
            double equivalentPoissonRatio, double thickness)
        {
            MaterialUtilities.CheckYoungModulus(youngModulus);
            MaterialUtilities.CheckPoissonRatio(poissonRatio);
            MaterialUtilities.CheckThickness(thickness);

            HomogeneousYoungModulus = youngModulus;
            HomogeneousEquivalentYoungModulus = equivalentYoungModulus;
            HomogeneousPoissonRatio = poissonRatio;
            HomogeneousEquivalentPoissonRatio = equivalentPoissonRatio;
            HomogeneousThickness = thickness;
        }

        public double GetYoungModulusAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        { return HomogeneousYoungModulus; }

        public double GetEquivalentYoungModulusAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        { return HomogeneousEquivalentYoungModulus; }

        public double GetPoissonRatioAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        { return HomogeneousPoissonRatio; }

        public double GetEquivalentPoissonRatioAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        { return HomogeneousEquivalentPoissonRatio; }

        public double GetThicknessAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        { return HomogeneousThickness; }
        
        public Matrix CalculateConstitutiveMatrixAt(INaturalPoint2D point, EvaluatedInterpolation2D interpolation)
        {
            var matrix = Matrix.CreateZero(3, 3);
            double eqE = HomogeneousEquivalentYoungModulus;
            double eqV = HomogeneousEquivalentPoissonRatio;
            double scalar = eqE / (1 - eqV * eqV);
            matrix[0, 0] = scalar;
            matrix[0, 1] = scalar * eqV;
            matrix[1, 0] = scalar * eqV;
            matrix[1, 1] = scalar;
            matrix[2, 2] = 0.5 * eqE / (1 + eqV);
            return matrix;
        }

        public void UpdateDistributions() { } // Do nothing for homogeneous preperties
        public void UpdateStrains() { } // Do nothing for elastic properties.
    }
}
