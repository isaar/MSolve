using ISAAR.MSolve.FEM.Interpolation;
using System;
using System.Collections.Generic;
using System.Text;

// TODO: move the node orderings somewhere more central
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
                //{ InterpolationLine2.UniqueInstance, 3 },  // 0 ---- 1
                //{ InterpolationLine2.UniqueInstance, 21 }, // 0 -- 2 -- 1

                                                             // 3 ---- 2
                                                             // |      |
                                                             // |      |
                { InterpolationQuad4.UniqueInstance, 9 },    // 0 ---- 1

                //                                             // 3 -- 6 -- 2
                //                                             // |         |
                //                                             // 7         5
                //                                             // |         |
                //{ InterpolationQuad8.UniqueInstance, 23 }    // 0 -- 4 -- 1

                                                             // 2
                                                             // | \
                                                             // |   \
                { InterpolationTri3.UniqueInstance, 5 }      // 0 --- 1

                //,
                //                                             // 2
                //                                             // | \
                //                                             // 5   4
                //                                             // |     \
                //{ InterpolationTri6.UniqueInstance, 22 }     // 0 - 3 - 1
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
