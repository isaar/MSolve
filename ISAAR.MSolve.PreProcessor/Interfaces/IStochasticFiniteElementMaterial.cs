using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IStochasticFiniteElementMaterial : IFiniteElementMaterial3D
    {
        IStochasticMaterialCoefficientsProvider CoefficientsProvider { get; set; }
        IMatrix2D<double> GetConstitutiveMatrix(double[] coordinates);
    }
}
