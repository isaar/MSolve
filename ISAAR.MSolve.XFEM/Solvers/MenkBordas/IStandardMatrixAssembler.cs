using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    interface IStandardMatrixAssembler
    {
        void BuildStandardMatrices(Model2D model, XClusterDofOrderer globalDofOrderer);
        DokRowMajor Ksc { get; }
    }
}
