using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;


namespace ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses
{
    public static class PrintUtilities
    {
        public static void WriteToFile(double[,] array, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                for (int j = 0; j < array.GetLength(1); ++j)
                {
                    writer.Write(array[i, j]);
                    writer.Write(' ');
                }
                writer.WriteLine();
            }
            writer.Flush();
        }

        public static void WriteToFile(SkylineMatrix2D Mat, int i1, int j1, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < 40; ++i)
            {
                for (int j = 0; j < 40; ++j)
                {
                    writer.Write(Mat[i + i1, j + j1]);
                    writer.Write(' ');
                }
                writer.WriteLine();
            }
            writer.Flush();
        }

        public static double[] ReadVector(string path)
        {

            var reader = new StreamReader(path);
            var lines = File.ReadLines(path).Count();
            double[] data = new double[lines];
            for (int i = 0; i < lines; ++i)
            {
                data[i] = Convert.ToDouble(reader.ReadLine());

            }
            reader.Close();
            return data;
        }

        public static int[] ReadIntVector(string path)
        {

            var reader = new StreamReader(path);
            var lines = File.ReadLines(path).Count();
            int[] data = new int[lines];
            for (int i = 0; i < lines; ++i)
            {
                data[i] = Convert.ToInt32(reader.ReadLine());

            }
            reader.Close();
            return data;
        }

        public static void WriteToFileVector(double[] array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }
            writer2.Flush();

        }

        public static void ConvertAndWriteToFileVector(double[][] array, string path3)
        {
            int length1 = array.GetLength(0);
            int length2 = array[0].GetLength(0);
            int vec_length = length1 * length2;
            double[] vector = new double[vec_length];
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                for (int j = 0; j < array[0].GetLength(0); ++j)
                {
                    vector[length2 * i + j] = array[i][j];
                }
            }
            WriteToFileVector(vector, path3);
        }

        public static void SeparateAndWriteToFile(double[,] array, string path_A, string path_B)
        {

            int length1 = array.GetLength(0);
            int length2 = array.GetLength(1);
            double[,] array_A;
            double[,] array_B;
            array_A = new double[40, length2];
            array_B = new double[length1 - 40, length2];
            // opou 24 --> length1-40
            for (int n = 0; n < 40; n++)
            {
                for (int p = 0; p < length2; p++)
                {
                    array_A[n, p] = array[n, p];
                }
            }
            for (int n = 0; n < length1 - 40; n++)
            {
                for (int p = 0; p < length2; p++)
                {
                    array_B[n, p] = array[n + 40, p];
                }
            }

            WriteToFile(array_A, path_A);
            WriteToFile(array_B, path_B);

        }
    }
}
