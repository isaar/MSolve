using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Materials
{
    public class StochasticElasticMaterial : IStochasticFiniteElementMaterial
    {
        public IStochasticMaterialCoefficientsProvider CoefficientsProvider { get; set; }

        private double _youngModulus;
        public double YoungModulus
        {
            get {return CoefficientsProvider.GetCoefficient(0, null); }
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
            return new double[] { stochasticYoungModulus,  PoissonRatio};
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

        public double[] Stresses => throw new NotImplementedException();

        public IMatrix2D ConstitutiveMatrix => throw new NotImplementedException();

        public void ResetModified()
        {
            throw new NotImplementedException();
        }

        public object Clone()
        {
            return new StochasticElasticMaterial (CoefficientsProvider) { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }
        #endregion

        #region IStochasticFiniteElementMaterial

        public IMatrix2D GetConstitutiveMatrix(double[] coordinates)
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



