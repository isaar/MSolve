using System.IO;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Input
{
    /// <summary>
    /// Tests for <see cref="Array1DReader"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class Array1DReaderTests
    {
        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);
        private static readonly double[] expectedArray = { 1.0, 2.0, 3.0, 4.0, 5.0 };

        [Fact]
        private static void TestArray1()
        {
            string path = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName 
                + @"\Resources\Array1DSingleLineAndLength.txt";
            var reader = new Array1DReader(true);
            comparer.AssertEqual(expectedArray, reader.ReadFile(path));
        }

        [Fact]
        private static void TestArray2()
        {
            string path = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\Array1DSingleLineAndLength.txt";
            var reader = new Array1DReader(true, ' ');
            comparer.AssertEqual(expectedArray, reader.ReadFile(path));
        }

        [Fact]
        private static void TestArray3()
        {
            string path = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\Array1DSingleLine.txt";
            var reader = new Array1DReader(false);
            comparer.AssertEqual(expectedArray, reader.ReadFile(path));
        }

        [Fact]
        private static void TestArray4()
        {
            string path = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\Array1DLength.txt";
            var reader = new Array1DReader(true);
            comparer.AssertEqual(expectedArray, reader.ReadFile(path));
        }

        [Fact]
        private static void TestArray5()
        {
            string path = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
                + @"\Resources\Array1DLength.txt";
            var reader = new Array1DReader(true, '\n');
            comparer.AssertEqual(expectedArray, reader.ReadFile(path));
        }
    }
}
