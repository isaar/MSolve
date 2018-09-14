using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Input
{
    /// <summary>
    /// Reads the index and value arrays of a skyline matrix from separate files.
    /// Authors: George Stavroulakis
    /// </summary>
    public class SkylineMatrixReader
    {
        public static SkylineMatrix ReadFromFiles(string valuesPath, string diagonalOffsetsPath)
        {
            FileStream fs = File.OpenRead(diagonalOffsetsPath);
            BinaryReader bw = new BinaryReader(fs);
            int length = 0;
            length = bw.ReadInt32();
            int[] diagOffsets = new int[length];
            for (int i = 0; i < diagOffsets.Length; i++) diagOffsets[i] = bw.ReadInt32();
            bw.Close();
            fs.Close();

            fs = File.OpenRead(valuesPath);
            bw = new BinaryReader(fs);
            length = bw.ReadInt32();
            double[] values = new double[length];
            for (int i = 0; i < values.Length; i++) values[i] = bw.ReadDouble();
            bw.Close();
            fs.Close();

            //string[] lines = File.ReadAllLines(path + "\\" + nameOnly + "-Ix" + ext);
            //rowIndex = new int[lines.Length];
            //for (int i = 0; i < lines.Length; i++) rowIndex[i] = Int32.Parse(lines[i]);

            //lines = File.ReadAllLines(path + "\\" + nameOnly + "-Data" + ext);
            //data = new T[lines.Length];
            //double[] mData = data as double[];
            //for (int i = 0; i < lines.Length; i++) mData[i] = Convert.ToDouble(lines[i]);

            return SkylineMatrix.CreateFromArrays(diagOffsets.Length - 1, values, diagOffsets, true, false);
        }
    }
}
