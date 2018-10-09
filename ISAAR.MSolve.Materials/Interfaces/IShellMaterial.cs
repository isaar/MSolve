using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IShellMaterial : IFiniteElementMaterial
    {
        double[] Stresses { get; }
        IMatrix2D ConstitutiveMatrix { get; }
        new IShellMaterial Clone();
        double[] NormalVectorV3 { get; set; } 
        double[] TangentVectorV1 { get; set; }
        void UpdateMaterial(double[] GLvec);
    }
}
