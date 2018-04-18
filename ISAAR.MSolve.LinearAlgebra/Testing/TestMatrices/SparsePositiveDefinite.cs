using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    public class SparsePositiveDefinite
    {
        public static readonly int order = 10;

        public static readonly double[,] matrix = {
            { 21.0,  1.0,  0.0,  4.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  1.0, 22.0,  2.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0 },
            {  0.0,  2.0, 23.0,  1.0,  3.0,  1.0,  0.0,  1.0,  0.0,  0.0 },
            {  4.0,  0.0,  1.0, 24.0,  2.0,  4.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  3.0,  2.0, 25.0,  5.0,  2.0,  0.0,  0.0,  1.0 },
            {  0.0,  0.0,  1.0,  4.0,  5.0, 26.0,  0.0,  0.0,  2.0,  3.0 },
            {  0.0,  1.0,  0.0,  0.0,  2.0,  0.0, 27.0,  3.0,  0.0,  0.0 },
            {  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  3.0, 28.0,  4.0,  2.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  2.0,  0.0,  4.0, 29.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  1.0,  3.0,  0.0,  2.0,  0.0, 30.0 }};

        public static readonly double[] skylineValues = {
            21.0, 22.0, 1.0,
            23.0, 2.0,
            24.0, 1.0, 0.0, 4.0,
            25.0, 2.0, 3.0,
            26.0, 5.0, 4.0, 1.0,
            27.0, 0.0, 2.0, 0.0, 0.0, 1.0,
            28.0, 3.0, 0.0, 0.0, 0.0, 1.0,
            29.0, 4.0, 0.0, 2.0,
            30.0, 0.0, 2.0, 0.0, 3.0, 1.0 };

        public static readonly int[] skylineDiagOffsets = { 0, 1, 3, 5, 9, 12, 16, 22, 28, 32, 38 };

        /// <summary>
        /// A = transpose(U) * U
        /// </summary>
        public static readonly double[,] choleskyU = {
            { 4.582575694955840, 0.218217890235992, 0.000000000000000,  0.872871560943970, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 4.685336802448780, 0.426863656622231, -0.040653681583070, 0.000000000000000, 0.000000000000000,  0.213431828311116,  0.000000000000000, 0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 4.776796773849092,  0.212978200107921, 0.628035929102890, 0.209345309700963, -0.019072674636530,  0.209345309700963, 0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  4.815712076375390, 0.387531897384781, 0.821356000941806,  0.002645268924128, -0.009258441234449, 0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 4.945239114569206, 0.920121731028097,  0.406644279935036, -0.025860920335721, 0.000000000000000,  0.202214691106420 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 4.943169515716910, -0.075324580738928, -0.002513728813410, 0.404598708104377,  0.569257813116235 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000,  5.175233591581246,  0.582455663587526, 0.005888860380148, -0.007603587481719 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  5.255107907216917, 0.760705036305723,  0.382692269290029 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 5.315787152855426, -0.098083711953027 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000,  5.429449618408797 }};

        public static readonly double[] lhs = { 7.5025, 9.8100, 2.3352, 0.9623, 3.8458, 5.0027, 5.7026, 9.7663, 4.9286, 4.0088 };
        public static readonly double[] rhs = {
            171.2117, 233.6955, 100.5983, 83.1428, 145.5027, 177.3672, 200.7707, 320.6314, 192.0000, 158.6505 };

        public static void CheckEquals()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var full = Matrix.CreateFromArray(matrix);
            var skyline = SkylineMatrix.CreateFromArrays(order, skylineValues, skylineDiagOffsets, true, true);
            if (skyline.Equals(full)) Console.WriteLine("Skyline.Equals() works fine.");
            else Console.WriteLine("Error in Skyline.Equals().");
            if (full.Equals(skyline)) Console.WriteLine("Matrix.Equals() works fine.");
            else Console.WriteLine("Error in Matrix.Equals().");
        }

        public static void CheckIndexing()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var skyline = SkylineMatrix.CreateFromArrays(order, skylineValues, skylineDiagOffsets, true, true);
            comparer.CheckMatrixEquality(matrix, skyline.CopyToFullMatrix().CopyToArray2D());
        }

        public static void CheckFactorization()
        {
            // Either my reconstruction of D*U or matlab's cholesky has limited precision. Probably mine, but I cannot figure why.
            //var comparer = new Comparer(Comparer.PrintMode.Always, 1e-1); 
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var skyline = SkylineMatrix.CreateFromArrays(order, skylineValues, skylineDiagOffsets, true, true);

            try
            {
                var factor = skyline.FactorCholesky(false);
                double[,] factorizedU = factor.GetFactorU().CopyToArray2D();
                comparer.CheckFactorizationCholesky(matrix, choleskyU, factorizedU);
            }
            catch (IndefiniteMatrixException)
            {
                var printer = new Printer();
                printer.PrintIndefiniteMatrix(matrix);
            }
        }

        public static void CheckMatrixVectorMult()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var skyline = SkylineMatrix.CreateFromArrays(order, skylineValues, skylineDiagOffsets, true, true);
            var x = Vector.CreateFromArray(lhs);
            Vector b = skyline.MultiplyRight(x, false);
            comparer.CheckMatrixVectorMult(matrix, lhs, rhs, b.CopyToArray());
        }

        public static void CheckSystemSolution()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var b = Vector.CreateFromArray(rhs);
            var skyline = SkylineMatrix.CreateFromArrays(order, skylineValues, skylineDiagOffsets, true, true);
            try
            {
                var factor = skyline.FactorCholesky(false);
                var x = factor.SolveLinearSystem(b);
                comparer.CheckSystemSolution(matrix, rhs, lhs, x.InternalData);
            }
            catch (IndefiniteMatrixException)
            {
                var printer = new Printer();
                printer.PrintIndefiniteMatrix(matrix);
            }
        }

        public static void Print()
        {
            var skyline = SkylineMatrix.CreateFromArrays(order, skylineValues, skylineDiagOffsets, true, true);
            INumericFormat storedDefault = FullMatrixWriter.NumericFormat;
            FullMatrixWriter.NumericFormat = new FixedPointFormat { MaxIntegerDigits = 2 };
            Console.WriteLine("Skyline (full) = ");
            (new FullMatrixWriter(skyline)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("Skyline (arrays) = ");
            (new RawArraysWriter(skyline, false)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("Skyline (sparse entries) = ");
            (new CoordinateTextFileWriter(skyline)).WriteToConsole();
            Console.WriteLine();
            FullMatrixWriter.NumericFormat = storedDefault;
        }
    }
}
