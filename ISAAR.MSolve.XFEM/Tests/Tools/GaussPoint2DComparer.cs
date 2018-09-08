using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    /// <summary>
    /// eta major, i.e. points with the same eta will be consecutive 
    /// </summary>
    class GaussPoint2DComparer : IComparer<GaussPoint2D>
    {
        private readonly ValueComparer valueComparer;

        public GaussPoint2DComparer(double tolerance)
        {
            this.valueComparer = new ValueComparer(tolerance);
        }

        public int Compare(GaussPoint2D point1, GaussPoint2D point2)
        {
            if (valueComparer.AreEqual(point1.Eta, point2.Eta))
            {
                if (valueComparer.AreEqual(point1.Xi, point2.Xi)) return 0;
                else if (point1.Xi < point2.Xi) return -1;
                else return 1;
            }
            else if (point1.Eta < point2.Eta) return -1;
            else return 1;
        }
    }
}
