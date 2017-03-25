using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.XFEM.Materials
{
    class ElasticMaterial2DPlainStrain : IFiniteElementMaterial2D
    {
        public double YoungModulus { get; }
        public double PoissonRatio { get; }
        public double Thickness { get; }

        public static ElasticMaterial2DPlainStrain Create(double youngModulus, double poissonRatio, double thickness)
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
            return new ElasticMaterial2DPlainStrain(youngModulus, poissonRatio, thickness);
        }

        private ElasticMaterial2DPlainStrain(double youngModulus, double poissonRatio, double thickness)
        {
            this.YoungModulus = youngModulus;
            this.PoissonRatio = poissonRatio;
            this.Thickness = thickness;
        }

        public Matrix2D CalculateConstitutiveMatrix() //This needs to be stored and passed as a MatrixView
        {
            double plainE = YoungModulus / (1 - PoissonRatio * PoissonRatio);
            double plainV = PoissonRatio / (1.0 - PoissonRatio); 
            double scalar = plainE / (1 - plainV * plainV);

            var matrix = new Matrix2D(3, 3);
            matrix[0, 0] = scalar;
            matrix[0, 1] = scalar * plainV;
            matrix[1, 0] = scalar * plainV;
            matrix[1, 1] = scalar;
            matrix[2, 2] = 0.5 * plainE / (1 + plainV);
            return matrix;
        }

        public IFiniteElementMaterial2D Clone()
        {
            return this; // For now it is immutable otherwise:
            //return new ElasticMaterial2DPlainStrain(YoungModulus, PoissonRatio, Thickness);
        }
    }
}
