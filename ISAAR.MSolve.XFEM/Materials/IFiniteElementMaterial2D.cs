using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.XFEM.Materials
{
    public interface IFiniteElementMaterial2D
    {
        double YoungModulus { get; }
        double PoissonRatio { get; }
        double Thickness { get; }
        Matrix2D<double> CalculateConstitutiveMatrix();
        IFiniteElementMaterial2D Clone();

        #region non linearity
        //bool Modified { get; }
        //void ResetModified();
        //void UpdateMaterial(double[] strains);
        //void ClearState();
        //void SaveState();
        //void ClearStresses();
        #endregion
    }
}
