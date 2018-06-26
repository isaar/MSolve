using ISAAR.MSolve.FEM.Interpolation;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging.VTK
{
    /// <summary>
    /// Cell used to represent VTK grids.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkCell2D
    {
        public static readonly IReadOnlyDictionary<IIsoparametricInterpolation2D, int> cellTypeCodes = 
            new Dictionary<IIsoparametricInterpolation2D, int>
            {
                { InterpolationQuad4.UniqueInstance, 9 }
            };

        public VtkCell2D(int code, IReadOnlyList<VtkPoint2D> vertices)
        {
            this.Code = code;
            this.Vertices = vertices;
        }

        public int Code { get; }
        public IReadOnlyList<VtkPoint2D> Vertices { get; }
    }
}
