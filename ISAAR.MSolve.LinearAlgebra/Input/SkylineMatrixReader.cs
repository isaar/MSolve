using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: this doesn't work for text files.
namespace ISAAR.MSolve.LinearAlgebra.Input
{
    /// <summary>
    /// Reads the index and value arrays of a skyline matrix from separate files or a single one.
    /// Authors: George Stavroulakis
    /// </summary>
    public static class SkylineMatrixReader
    {
        /// <summary>
        /// Reads the Skyline values and diagonal offsets arrays from 2 different files. The first entry in each file must be 
        /// the length of the corresponding array.
        /// </summary>
        /// <param name="valuesPath">
        /// The absolute path of an array containing the values array of the Skyline format. The first entry must be its entry, 
        /// which is equal to the number of nonzero entries.
        /// </param>
        /// <param name="diagonalOffsetsPath">
        /// The absolute path of an array containing the diagonal offsets array of the Skyline format. The first entry must be 
        /// its entry, which is equal to the number of rows/columns + 1.
        /// </param>
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

            return SkylineMatrix.CreateFromArrays(diagOffsets.Length - 1, values, diagOffsets, true, false);
        }

        public static SkylineMatrix ReadFromSimilarlyNamedFiles(string commonNamePart)
        {
            string path = Path.GetDirectoryName(commonNamePart);
            string nameOnly = Path.GetFileNameWithoutExtension(commonNamePart);
            string ext = Path.GetExtension(commonNamePart);

            FileStream fs = File.OpenRead(path + "\\" + nameOnly + "-Ix" + ext);
            BinaryReader bw = new BinaryReader(fs);
            int length = 0;
            length = bw.ReadInt32();
            var diagOffsets = new int[length];
            for (int i = 0; i < diagOffsets.Length; i++) diagOffsets[i] = bw.ReadInt32();
            bw.Close();
            fs.Close();

            fs = File.OpenRead(path + "\\" + nameOnly + "-Data" + ext);
            bw = new BinaryReader(fs);
            length = bw.ReadInt32();
            var values = new double[length];
            for (int i = 0; i < values.Length; i++) values[i] = bw.ReadDouble();
            bw.Close();
            fs.Close();

            return SkylineMatrix.CreateFromArrays(diagOffsets.Length - 1, values, diagOffsets, true, false);

            //string[] lines = File.ReadAllLines(path + "\\" + nameOnly + "-Ix" + ext);
            //rowIndex = new int[lines.Length];
            //for (int i = 0; i < lines.Length; i++) rowIndex[i] = Int32.Parse(lines[i]);

            //lines = File.ReadAllLines(path + "\\" + nameOnly + "-Data" + ext);
            //data = new T[lines.Length];
            //double[] mData = data as double[];
            //for (int i = 0; i < lines.Length; i++) mData[i] = Convert.ToDouble(lines[i]);
        }
    }
}
