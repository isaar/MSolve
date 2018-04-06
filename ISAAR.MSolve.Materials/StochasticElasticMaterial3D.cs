using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;


namespace ISAAR.MSolve.FEM.Materials
{
    public class StochasticElasticMaterial3D : IStochasticIsotropicContinuumMaterial3D, IIsotropicContinuumMaterial3D//why do we need both interfaces?? IIsotropicFiniteElementMaterial3D or IStochasticIsotropicFiniteElementMaterial will be fine. 
    {
        private IStochasticMaterialCoefficientsProvider coefficientsProvider;
        private readonly double[] strains = new double[6];
        private readonly double[] stresses = new double[6];
        //private double[,] constitutiveMatrix = null;
        public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public double[] Coordinates { get; set; }

        public StochasticElasticMaterial3D(IStochasticMaterialCoefficientsProvider coefficientsProvider)
        {
            this.coefficientsProvider = coefficientsProvider;
        }

        private double[,] GetConstitutiveMatrixInternal(double[] coordinates)
        {
            //double variation = coefficientsProvider.GetCoefficient(YoungModulus, coordinates);
            //double fE1 = variation * YoungModulus / (1.0 + PoissonRatio);
            double stochasticYoungModulus = coefficientsProvider.GetCoefficient(YoungModulus, coordinates);
            double fE1 = stochasticYoungModulus / (1.0 + PoissonRatio);
            double fE2 = fE1 * PoissonRatio / (1.0 - 2.0 * PoissonRatio);
            double fE3 = fE1 + fE2;
            double fE4 = fE1 * 0.5;
            double[,] afE = new double[6, 6];
            afE[0, 0] = fE3;
            afE[0, 1] = fE2;
            afE[0, 2] = fE2;
            afE[1, 0] = fE2;
            afE[1, 1] = fE3;
            afE[1, 2] = fE2;
            afE[2, 0] = fE2;
            afE[2, 1] = fE2;
            afE[2, 2] = fE3;
            afE[3, 3] = fE4;
            afE[4, 4] = fE4;
            afE[5, 5] = fE4;

            //Vector<double> s = (new Matrix2D<double>(afE)) * (new Vector<double>(strains));
            //s.Data.CopyTo(stresses, 0);

            return afE;
        }

        #region IFiniteElementMaterial Members

        public int ID
        {
            get { return 1; }
        }

        public bool Modified
        {
            get { return false; }
        }

        public void ResetModified()
        {
        }

        #endregion

        #region IFiniteElementMaterial3D Members

        public StressStrainVectorContinuum3D Stresses { get { return new StressStrainVectorContinuum3D(new double[6]); } }

        public ElasticityTensorContinuum3D ConstitutiveMatrix
        {
            get
            {
                throw new NotImplementedException();
                //if (constitutiveMatrix == null) UpdateMaterial(new double[6]);
                //return new Matrix2D<double>(constitutiveMatrix);
            }
        }

        public void UpdateMaterial(StressStrainVectorContinuum3D strains)
        {
            throw new NotImplementedException();

            //strains.CopyTo(this.strains, 0);
            //constitutiveMatrix = GetConstitutiveMatrix();
        }

        public void ClearState()
        {
        }

        public void SaveState()
        {
            throw new NotImplementedException();
        }

        public void ClearStresses()
        {
            throw new NotImplementedException();
        }

        #endregion

        #region ICloneable Members

        public object Clone()
        {
            return new StochasticElasticMaterial3D(coefficientsProvider) { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

        #endregion

        #region IStochasticFiniteElementMaterial Members
        public IStochasticMaterialCoefficientsProvider CoefficientsProvider
        {
            get { return coefficientsProvider; }
            set { coefficientsProvider = value; }
        }

        public ElasticityTensorContinuum3D GetConstitutiveMatrix(double[] coordinates)
        {
            return new ElasticityTensorContinuum3D(GetConstitutiveMatrixInternal(coordinates));
        }

        #endregion
    }
}
