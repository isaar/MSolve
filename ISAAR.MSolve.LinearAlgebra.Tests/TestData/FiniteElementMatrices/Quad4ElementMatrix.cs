using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData.FiniteElementMatrices
{
    /// <summary>
    /// Stiffness matrix of a 2D qudrilateral element with 4 nodes, taken from Abaqus. The Quad4 element is arbitrary 
    /// (not a rectangle).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    class Quad4ElementMatrix
    {
        internal static double[,] UnconstrainedStiffness => new double[,]
        {
            { 181603.19122884000, -89089.52288016200, -4991.10465483030, 7519.91442892909, -126789.16986973000, 70773.21914644800, -49822.91670427800, 10796.38930478500 },
            { -89089.52288016200, 210532.99196037000, 13289.14519816000, -85869.32985046300, 70773.21914644800, -146868.68100448000, 5027.15853555400, 22205.01889457400 },
            { -4991.10465483030, 13289.14519816000, 73155.65851293700, 4056.21872780570, -66433.22620703800, 10577.55360637600, -1731.32765106900, -27922.91753234200 },
            { 7519.91442892909, -85869.32985046300, 4056.21872780570, 85360.61861633800, 16346.78437560700, 2912.80345339600, -27922.91753234200, -2404.09221927070 },
            { -126789.16986973000, 70773.21914644800, -66433.22620703800, 16346.78437560700, 200705.04715702000, -95472.47721160800, -7482.65108024430, 8352.47368955240 },
            { 70773.21914644800, -146868.68100448000, 10577.55360637600, 2912.80345339600, -95472.47721160800, 232719.03971773000, 14121.70445878300, -88763.16216664000 },
            { -49822.91670427800, 5027.15853555400, -1731.32765106900, -27922.91753234200, -7482.65108024430, 14121.70445878300, 59036.89543559100, 8774.05453800470 },
            { 10796.38930478500, 22205.01889457400, -27922.91753234200, -2404.09221927070, 8352.47368955240, -88763.16216664000, 8774.05453800470, 68962.23549133600 }
        }; // from Abaqus

        internal static double[,] ClampedStiffness
        {
            get
            {
                var unconstrained = Matrix.CreateFromArray(UnconstrainedStiffness);

                // constrain x, y dofs of node 0 and y dof of node 1
                var freeDofs = new int[] { 2, 4, 5, 6, 7 };
                return unconstrained.GetSubmatrix(freeDofs, freeDofs).CopyToArray2D();
            }
        }

        internal static double[,] PinnedStiffness
        {
            get
            {
                var unconstrained = Matrix.CreateFromArray(UnconstrainedStiffness);
                int order = unconstrained.NumColumns;

                // constrain x, y dofs of node 0
                return unconstrained.GetSubmatrix(2, order, 2, order).CopyToArray2D();
            }
        }
    }
}
