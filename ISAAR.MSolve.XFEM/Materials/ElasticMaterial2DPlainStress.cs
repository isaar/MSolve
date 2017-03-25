using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.XFEM.Materials
{
    public class ElasticMaterial2DPlainStress : IFiniteElementMaterial2D
    {
        public double YoungModulus { get; }
        public double PoissonRatio { get; }
        public double Thickness { get; }

        public static ElasticMaterial2DPlainStress Create(double youngModulus, double poissonRatio, double thickness)
        {
            if (youngModulus <= 0.0)
            {
                throw new ArgumentException("Young's modulus must be positive but was: " + youngModulus);
            }
            if (poissonRatio < -1.0 || poissonRatio > 0.5)
            {
                throw new ArgumentException("Poisson's ratio must be in the range [-1, 0.5] but was: " + poissonRatio);
            }
            if (thickness <= 0.0)
            {
                throw new ArgumentException("Thickness must be positive but was: " + thickness);
            }
            return new ElasticMaterial2DPlainStress(youngModulus, poissonRatio, thickness);
        }

        private ElasticMaterial2DPlainStress(double youngModulus, double poissonRatio, double thickness)
        {
            this.YoungModulus = youngModulus;
            this.PoissonRatio = poissonRatio;
            this.Thickness = thickness;
        }
        
        public Matrix2D CalculateConstitutiveMatrix() //This needs to be stored and passed as a MatrixView
        {
            var matrix = new Matrix2D(3, 3);
            double scalar = YoungModulus / (1 - PoissonRatio * PoissonRatio);
            matrix[0, 0] = scalar;
            matrix[0, 1] = scalar * PoissonRatio;
            matrix[1, 0] = scalar * PoissonRatio;
            matrix[1, 1] = scalar;
            matrix[2, 2] = 0.5 * YoungModulus / (1 + PoissonRatio);
            return matrix;
        }

        public IFiniteElementMaterial2D Clone()
        {
            return this; // For now it is immutable otherwise:
            //return new ElasticMaterial2DPlainStress(YoungModulus, PoissonRatio, Thickness);
        }

    }
}
