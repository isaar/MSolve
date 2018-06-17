using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace ISAAR.MSolve.FEM.Materials
{
    public class StochasticElasticMaterial : IStochasticContinuumMaterial3D
    {
        public IStochasticMaterialCoefficientsProvider CoefficientsProvider { get; set; }

        private double _youngModulus;
        public double YoungModulus
        {
            get { return CoefficientsProvider.GetCoefficient(0, null); }
            set { _youngModulus = value; }
        }

        public double PoissonRatio { get; set; }
        public double[] Coordinates { get; set; }

        public StochasticElasticMaterial(IStochasticMaterialCoefficientsProvider coefficientsProvider)
        {
            this.CoefficientsProvider = coefficientsProvider;
        }

        public double[] GetStochasticMaterialProperties(double[] coordinates)
        {
            double stochasticYoungModulus = CoefficientsProvider.GetCoefficient(YoungModulus, coordinates);
            return new double[] { stochasticYoungModulus, PoissonRatio };
        }

        #region IFiniteElementMaterialMembers

        public int ID
        {
            get { return 1; }
        }

        public bool Modified
        {
            get { return false; }
        }

        public StressStrainVectorContinuum3D Stresses => throw new NotImplementedException();

        public ElasticityTensorContinuum3D ConstitutiveMatrix => throw new NotImplementedException();

        public void ResetModified()
        {
            throw new NotImplementedException();
        }

        public object Clone()
        {
            return new StochasticElasticMaterial(CoefficientsProvider) { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }
        #endregion

        #region IStochasticFiniteElementMaterial

        public ElasticityTensorContinuum3D GetConstitutiveMatrix(double[] coordinates)
        {
            throw new NotImplementedException();
        }

        public void UpdateMaterial(StressStrainVectorContinuum3D strains)
        {
            throw new NotImplementedException();
        }

        public void ClearState()
        {
            throw new NotImplementedException();
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

    }
}



