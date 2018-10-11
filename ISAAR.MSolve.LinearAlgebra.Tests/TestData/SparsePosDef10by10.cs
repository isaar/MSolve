using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// A 10-by-10 sparse matrix that is symmetric positive definite.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SparsePosDef10by10
    {
        internal const int order = 10;

        internal static readonly double[,] matrix = {
            { 21.0,  1.0,  0.0,  4.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  1.0, 22.0,  2.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0 },
            {  0.0,  2.0, 23.0,  1.0,  3.0,  1.0,  0.0,  1.0,  0.0,  0.0 },
            {  4.0,  0.0,  1.0, 24.0,  2.0,  4.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  3.0,  2.0, 25.0,  5.0,  2.0,  0.0,  0.0,  1.0 },
            {  0.0,  0.0,  1.0,  4.0,  5.0, 26.0,  0.0,  0.0,  2.0,  3.0 },
            {  0.0,  1.0,  0.0,  0.0,  2.0,  0.0, 27.0,  3.0,  0.0,  0.0 },
            {  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  3.0, 28.0,  4.0,  2.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  2.0,  0.0,  4.0, 29.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  1.0,  3.0,  0.0,  2.0,  0.0, 30.0 }
        };

        internal static readonly double[] skylineValues = {
            21.0, 22.0, 1.0,
            23.0, 2.0,
            24.0, 1.0, 0.0, 4.0,
            25.0, 2.0, 3.0,
            26.0, 5.0, 4.0, 1.0,
            27.0, 0.0, 2.0, 0.0, 0.0, 1.0,
            28.0, 3.0, 0.0, 0.0, 0.0, 1.0,
            29.0, 4.0, 0.0, 2.0,
            30.0, 0.0, 2.0, 0.0, 3.0, 1.0
        };

        internal static readonly int[] skylineDiagOffsets = { 0, 1, 3, 5, 9, 12, 16, 22, 28, 32, 38 };

        /// <summary>
        /// A = transpose(U) * U
        /// </summary>
        internal static readonly double[,] choleskyU = {
            { 4.582575694955840, 0.218217890235992, 0.000000000000000,  0.872871560943970, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 4.685336802448780, 0.426863656622231, -0.040653681583070, 0.000000000000000, 0.000000000000000,  0.213431828311116,  0.000000000000000, 0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 4.776796773849092,  0.212978200107921, 0.628035929102890, 0.209345309700963, -0.019072674636530,  0.209345309700963, 0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  4.815712076375390, 0.387531897384781, 0.821356000941806,  0.002645268924128, -0.009258441234449, 0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 4.945239114569206, 0.920121731028097,  0.406644279935036, -0.025860920335721, 0.000000000000000,  0.202214691106420 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 4.943169515716910, -0.075324580738928, -0.002513728813410, 0.404598708104377,  0.569257813116235 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000,  5.175233591581246,  0.582455663587526, 0.005888860380148, -0.007603587481719 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  5.255107907216917, 0.760705036305723,  0.382692269290029 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 5.315787152855426, -0.098083711953027 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000,  5.429449618408797 }
        };

        internal static readonly double[] lhs = {
            7.5025, 9.8100, 2.3352, 0.9623, 3.8458, 5.0027, 5.7026, 9.7663, 4.9286, 4.0088 };

        internal static readonly double[] rhs = {
            171.2117, 233.6955, 100.5983, 83.1428, 145.5027, 177.3672, 200.7707, 320.6314, 192.0000, 158.6505 };

        internal static readonly int[] matlabPermutationAMD = { 0, 1, 8, 9, 7, 2, 3, 4, 5, 6 };

        /// <summary>
        /// Enforce indices order:
        /// First group: 0, 1, 7. Second group: 3, 5, 6, 9. Third group: 2, 4, 8.
        /// The new positions of these indices should be grouped, such as first group is before second before third.
        /// </summary>
        internal static readonly int[] constraintsCAMD = { 0, 0, 2, 1, 2, 1, 1, 0, 2, 1 };

        internal static readonly int[] permutationCAMD = { 0, 1, 7, 3, 6, 9, 5, 2, 4, 8 };

        internal static readonly string fullFormatPath = 
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName 
            + @"\Resources\SparsePosDef10by10_FullFormat.txt";

        internal static readonly string skylineArraysPath =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparsePosDef10by10_SkylineArrays.txt";

        internal static readonly string coordinateFormatPath =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparsePosDef10by10_CoordinateFormat.txt";

        internal static readonly string patternPath =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
            + @"\Resources\SparsePosDef10by10_Pattern.txt";
    }
}
