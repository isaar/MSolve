using System;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.FEM.Materials
{
    public class StochasticElasticMaterial_v2 : IStochasticContinuumMaterial3D_v2
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

        public StochasticElasticMaterial_v2(IStochasticMaterialCoefficientsProvider coefficientsProvider)
        {
            this.CoefficientsProvider = coefficientsProvider;
        }

        public double[] GetStochasticMaterialProperties(double[] coordinates)
        {
            double stochasticYoungModulus = CoefficientsProvider.GetCoefficient(YoungModulus, coordinates);
            return new double[] { stochasticYoungModulus, PoissonRatio };
        }

        #region IFiniteElementMaterialMembers

        public int ID => 1;

        public bool Modified => false;

        public double[] Stresses => throw new NotImplementedException();

        public IMatrixView ConstitutiveMatrix => throw new NotImplementedException();

        public void ResetModified()
        {
            throw new NotImplementedException();
        }

        public object Clone() => new StochasticElasticMaterial_v2(CoefficientsProvider)
        {
            YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio
        };
        #endregion

        #region IStochasticFiniteElementMaterial

        public IMatrixView GetConstitutiveMatrix(double[] coordinates)
        {
            throw new NotImplementedException();
        }

        public void UpdateMaterial(double[] strains)
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



