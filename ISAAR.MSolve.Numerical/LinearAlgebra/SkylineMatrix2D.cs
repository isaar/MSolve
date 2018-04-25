using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.IO;
using System.Globalization;
using ISAAR.MSolve.Numerical.Interfaces;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class SkylineMatrix2D : IMatrix2D, ICloneable, ILinearlyCombinable<SkylineMatrix2D>, ILinearlyCombinable, ISolveable, IFileReadable, IFileWriteable
    {
        private bool isFactorized;
        private double[] data;
        private int[] rowIndex;
        private double nullT = 0;

        public SkylineMatrix2D(int[] rowIndex)
        {
            this.rowIndex = rowIndex;
            data = rowIndex.Length > 0 ? new double[rowIndex[rowIndex.Length - 1]] : new double[0];
        }

        public SkylineMatrix2D(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            if (matrix.GetLength(1) != rows) throw new ArgumentException("Matrix must be square.");

            rowIndex = new int[rows + 1];
            for (int i = 0; i < rows; i++) rowIndex[i + 1] = rowIndex[i] + i + 1;
            data = new double[rowIndex[rowIndex.Length - 1]];

            int pos = 0;
            for (int j = 0; j < rows; j++)
                for (int i = j; i >= 0; i--)
                {
                    data[pos] = matrix[i, j];
                    pos++;
                }
        }

        public static SkylineMatrix2D Empty(int rows)
        {
            int[] rowIndex = new int[rows + 1];
            return new SkylineMatrix2D(rowIndex);
        }

        public int[] RowIndex
        {
            get { return rowIndex; }
        }

        public double[] Data
        {
            get { return data; }
        }

        public bool IsFactorized
        {
            get { return isFactorized; }
        }

        public void Factorize(double tolerance, List<IVector> zems, List<int> zemCols)
        {
            if (isFactorized)
                throw new InvalidOperationException("Matrix is already factorized.");
            double[] d = data;

            int kFix = 0;
            for (int n = 0; n < Rows; n++)
            {
                int KN = rowIndex[n];
                int KL = KN + 1;
                int KU = rowIndex[n + 1] - 1;
                int KH = KU - KL;
                if (KH < 0) continue;

                int K;
                if (KH > 0)
                {
                    K = n - KH;
                    //int IC = 0;
                    int KLT = KU;
                    for (int j = 0; j < KH; j++)
                    {
                        //IC++;
                        KLT--;
                        int KI = rowIndex[K];
                        int ND = rowIndex[K + 1] - KI - 1;
                        if (ND > 0)
                        {
                            //int KK = Math.Min(IC, ND);
                            int KK = Math.Min(j + 1, ND);
                            double C = 0;
                            for (int l = 1; l <= KK; l++)
                                C += d[KI + l] * d[KLT + l];
                            d[KLT] -= C;
                        }
                        K++;
                    }
                }
                K = n;
                double B = 0;
                for (int KK = KL; KK <= KU; KK++)
                {
                    K--;
                    int KI = rowIndex[K];
                    //if (d[KI] == 0) throw new InvalidOperationException(String.Format("Zero element in diagonal at index {0}.", KI));
                    if (Math.Abs(d[KI]) < tolerance) throw new InvalidOperationException(String.Format("Near-zero element in diagonal at index {0}.", KI));
                    double C = d[KK] / d[KI];
                    B += C * d[KK];
                    d[KK] = C;
                }
                d[KN] -= B;

                if (Math.Abs(d[KN]) < tolerance)
                {
                    d[KN] = 1;
                    zemCols.Add(n);
                    int j1 = n;
                    zems.Add(new Vector(Rows));
                    zems[kFix][j1] = 1;
                    for (int i1 = KN + 1; i1 <= KU; i1++)
                    {
                        j1--;
                        zems[kFix][j1] = -d[i1];
                        d[i1] = 0;
                    }
                    for (int irest = n + 1; irest < Rows; irest++)
                    {
                        int m1 = rowIndex[irest] + irest - n;
                        if (m1 <= rowIndex[irest + 1]) d[m1] = 0;
                    }
                    kFix++;
                }
            }
            isFactorized = true;
            if (Rows < 2) return;

            for (int ifl = 0; ifl < kFix; ifl++)
            {
                int n = Rows - 1;
                for (int l = 1; l < Rows; l++)
                {
                    int KL = rowIndex[n] + 1;
                    int KU = rowIndex[n + 1] - 1;
                    if (KU - KL >= 0)
                    {
                        int k = n;
                        for (int KK = KL; KK <= KU; KK++)
                        {
                            k--;
                            zems[ifl][k] -= d[KK] * zems[ifl][n];
                        }
                    }
                    n--;
                }
            }
        }

        public void Solve(IVector f, IVector result)
        {
            //var e = DateTime.Now;
            SkylineMatrix2D K = this;
            if (!K.isFactorized) throw new InvalidOperationException("Cannot solve if matrix is not factorized.");
            if (K.Rows != f.Length) throw new ArgumentException("Matrix and vector size mismatch.");
            double[] d = K.Data;
            //double[] result = new double[K.Rows];
            result.CopyFrom(0, f.Length, f, 0);
            // f.CopyTo(result, 0);

            // RHS vector reduction
            int n;
            for (n = 0; n < K.Rows; n++)
            {
                int KL = K.RowIndex[n] + 1;
                int KU = K.RowIndex[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    double C = 0;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        C += d[KK] * result[k];
                    }
                    result[n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < K.Rows; n++) result[n] /= d[K.RowIndex[n]];

            n = K.Rows - 1;
            for (int l = 1; l < K.Rows; l++)
            {
                int KL = K.RowIndex[n] + 1;
                int KU = K.RowIndex[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        result[k] -= d[KK] * result[n];
                    }
                }
                n--;
            }
            //var x = new List<TimeSpan>();
            //x.Add(DateTime.Now - e);
        }

        //public static double[] operator /(SkylineMatrix2D K, IVector f)
        //{
        //    if (!K.isFactorized) throw new InvalidOperationException("Cannot solve if matrix is not factorized.");
        //    if (K.Rows != f.Length) throw new ArgumentException("Matrix and vector size mismatch.");
        //    double[] result = new double[K.Rows];
        //    K.Solve(f, result);
        //    return result;
        //}

        public void Multiply(IVector vIn, double[] vOut, double scaleFactor, int vInStartIndex, int vOutStartIndex, bool clearvOut)
        {
            //if (isFactorized) throw new InvalidOperationException("Cannot multiply if matrix is factorized.");
            //if (!(typeof(T) == typeof(double))) throw new InvalidOperationException("Cannot multiply for types other than double");
            //if (Rows != vIn.Length) throw new InvalidOperationException("Matrix and vector size mismatch.");

            //Array.Clear(vOut, 0, vOut.Length);
            //double[] d = data as double[];
            //int pos = 0;
            //for (int i = 0; i < Rows; i++)
            //{
            //    int height = rowIndex[i + 1] - rowIndex[i];
            //    if (height <= 0) continue;
            //    vOut[i] += d[pos] * vIn[i];
            //    pos++;
            //    for (int j = 0; j < height - 1; j++)
            //    {
            //        int row = i - j - 1;
            //        vOut[row] += d[pos] * vIn[i];
            //        vOut[i] += d[pos] * vIn[row];
            //        pos++;
            //    }
            //}

            //            if (isFactorized) throw new InvalidOperationException("Cannot multiply if matrix is factorized.");
            //if (Rows < vIn.Length - vInStartIndex) throw new InvalidOperationException("Matrix and vector size mismatch.");

            if (clearvOut) Array.Clear(vOut, 0, vOut.Length);
            double[] d = data as double[];
            int pos = 0;
            for (int i = 0; i < Rows; i++)
            {
                int height = rowIndex[i + 1] - rowIndex[i];
                if (height <= 0) continue;

                vOut[i + vOutStartIndex] += scaleFactor * d[pos] * vIn[i + vInStartIndex];
                pos++;
                for (int j = 0; j < height - 1; j++)
                {
                    int row = i - j - 1;
                    vOut[row + vOutStartIndex] += scaleFactor * d[pos] * vIn[i + vInStartIndex];
                    vOut[i + vOutStartIndex] += scaleFactor * d[pos] * vIn[row + vInStartIndex];
                    pos++;
                }
            }

            //! Multiply SKYLINE K with vector pfVector and store to pfResult
            //iDataPos = 1
            //Do I = 1, iEqNo
            //    iColHeight = atMatrix.aiColumnIx(I + 1) - atMatrix.aiColumnIx(I)
            //    if (iColHeight.GT.0) then
            //        afResult(I) = afResult(I) + atMatrix.afData(iDataPos) * afVector(I);
            //        iDataPos = iDataPos + 1
            //        Do J = 1, iColHeight - 1
            //            iRow = I - J
            //            afResult(iRow) = afResult(iRow) + atMatrix.afData(iDataPos) * afVector(I)
            //            afResult(I) = afResult(I) + atMatrix.afData(iDataPos) * afVector(iRow)
            //            iDataPos = iDataPos + 1
            //        end Do
            //    end if
            //end Do
        }

        public void Multiply(IVector vIn, double[] vOut)
        {
            Multiply(vIn, vOut, 1.0, 0, 0, true);
        }

        public void Scale(double scale)
        {
            double[] mData = data;
            for (int i = 0; i < mData.Length; i++) mData[i] *= scale;
        }

        public static double[] operator *(SkylineMatrix2D K, IVector v)
        {
            if (K.isFactorized) throw new InvalidOperationException("Cannot multiply if matrix is factorized.");
            if (K.Rows != v.Length) throw new ArgumentException("Matrix and vector size mismatch.");
            double[] result = new double[K.Rows];
            K.Multiply(v, result);
            return result;
        }

        private void LinearCombinationInternal(IList<double> coefficients, IList<SkylineMatrix2D> matrices)
        {
            foreach (var matrix in matrices)
                if (matrix.Rows != this.Rows)
                    throw new ArgumentException("Matrices do not have the same size.");

            for (int i = 0; i < rowIndex.Length - 1; i++)
            {
                int currentMaxHeight = 0;
                foreach (var matrix in matrices)
                    currentMaxHeight = Math.Max(currentMaxHeight, matrix.RowIndex[i + 1] - matrix.RowIndex[i]);
                if (currentMaxHeight > rowIndex[i + 1] - rowIndex[i])
                    throw new InvalidOperationException("Current matrix does not have enough storage capacity for the requested linear combination.");
            }

            double temp;
            double[] d = data as double[];
            for (int i = 0; i < rowIndex.Length - 1; i++)
            {
                for (int j = rowIndex[i]; j < rowIndex[i + 1]; j++)
                {
                    temp = 0;
                    for (int k = 0; k < coefficients.Count; k++)
                    {
                        int pos = j - rowIndex[i];
                        if (pos < matrices[k].RowIndex[i + 1] - matrices[k].RowIndex[i])
                            temp += coefficients[k] * matrices[k].Data[matrices[k].RowIndex[i] + pos];
                    }
                    d[j] = temp;
                }
            }

            //for (int i = 0; i < data.Length; i++)
            //{
            //    temp = 0;
            //    for (int j = 0; j < coefficients.Count; j++) temp += coefficients[j] * matrices[j].Data[i];
            //    d[i] = temp;
            //}
        }

        #region IMatrix2D<T> Members

        public int Rows
        {
            get { return rowIndex.Length - 1; }
        }

        public int Columns
        {
            get { return rowIndex.Length - 1; }
        }

        public double this[int x, int y]
        {
            get
            {
                int r = x;
                int c = y;
                if (x < y)
                {
                    r = y;
                    c = x;
                }
                int row1Start = rowIndex[r];
                int row1Stop = rowIndex[r + 1];
                int cols1 = row1Stop - row1Start;
                int minCol1 = r - cols1 + 1;

                if (c >= minCol1)
                {
                    int offset = r - c;
                    int pos = row1Start + offset;
                    return data[pos];
                }

                return nullT;
                //if (x == y)
                //    return data[rowIndex[x]];

                //throw new NotImplementedException("Operator implemented only for diagonal elements.");
            }
            set
            {
                int r = x;
                int c = y;
                if (x < y)
                {
                    r = y;
                    c = x;
                }
                int row1Start = rowIndex[r];
                int row1Stop = rowIndex[r + 1];
                int cols1 = row1Stop - row1Start;
                int minCol1 = r - cols1 + 1;

                if (c >= minCol1)
                {
                    int offset = r - c;
                    int pos = row1Start + offset;
                    data[pos] = value;
                }
                else
                    throw new ArgumentException("Specified position is not indexed in the skyline format.");
            }
        }

        public void LinearCombination(IList<double> coefficients, IList<SkylineMatrix2D> matrices)
        {
            this.isFactorized = false;
            List<SkylineMatrix2D> m = new List<SkylineMatrix2D>(matrices.Count);
            foreach (var matrix in matrices) m.Add(matrix);
            LinearCombinationInternal((IList<double>)coefficients, m);
        }

        public void LinearCombination(IList<double> coefficients, IList<IMatrix2D> matrices)
        {
            this.isFactorized = false;
            List<SkylineMatrix2D> m = new List<SkylineMatrix2D>(matrices.Count);
            foreach (var matrix in matrices)
            {
                if (matrix is SkylineMatrix2D) m.Add((SkylineMatrix2D)matrix);
                else throw new InvalidCastException("Cannot linearly combine skyline matrix with other matrix types.");
            }
            LinearCombinationInternal((IList<double>)coefficients, m);
        }

        public void WriteToFile(string name)
        {
            double[] mData = data as double[];

            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            //// Binary
            //FileStream fs = File.Create(path + "\\" + nameOnly + "-Ix" + ext, 8192, FileOptions.SequentialScan);
            //BinaryWriter bw = new BinaryWriter(fs);
            //bw.Write(rowIndex.Length);
            //for (int i = 0; i < rowIndex.Length; i++) bw.Write(rowIndex[i]);
            //bw.Close();
            //fs.Close();

            //fs = File.Create(path + "\\" + nameOnly + "-Data" + ext, 8192, FileOptions.SequentialScan);
            //bw = new BinaryWriter(fs);
            //bw.Write(mData.Length);
            //for (int i = 0; i < mData.Length; i++) bw.Write(mData[i]);
            //bw.Close();
            //fs.Close();

            // ASCII
            StreamWriter sw = new StreamWriter(path + "\\" + nameOnly + "-Ix" + ext);
            foreach (int i in rowIndex) sw.WriteLine(i);
            sw.Close();

            sw = new StreamWriter(path + "\\" + nameOnly + "-Data" + ext);
            foreach (double d in mData)
                sw.WriteLine(d.ToString("g17", new CultureInfo("en-US", false).NumberFormat));
            sw.Close();
        }

        public void ReadFromFile(string name)
        {
            isFactorized = false;
            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            FileStream fs = File.OpenRead(path + "\\" + nameOnly + "-Ix" + ext);
            BinaryReader bw = new BinaryReader(fs);
            int length = 0;
            length = bw.ReadInt32();
            rowIndex = new int[length];
            for (int i = 0; i < rowIndex.Length; i++) rowIndex[i] = bw.ReadInt32();
            bw.Close();
            fs.Close();

            fs = File.OpenRead(path + "\\" + nameOnly + "-Data" + ext);
            bw = new BinaryReader(fs);
            length = bw.ReadInt32();
            data = new double[length];
            double[] mData = data as double[];
            for (int i = 0; i < mData.Length; i++) mData[i] = bw.ReadDouble();
            bw.Close();
            fs.Close();

            //string[] lines = File.ReadAllLines(path + "\\" + nameOnly + "-Ix" + ext);
            //rowIndex = new int[lines.Length];
            //for (int i = 0; i < lines.Length; i++) rowIndex[i] = Int32.Parse(lines[i]);

            //lines = File.ReadAllLines(path + "\\" + nameOnly + "-Data" + ext);
            //data = new T[lines.Length];
            //double[] mData = data as double[];
            //for (int i = 0; i < lines.Length; i++) mData[i] = Convert.ToDouble(lines[i]);
        }

        #endregion

        #region ICloneable Members

        public object Clone()
        {
            SkylineMatrix2D clone = new SkylineMatrix2D(this.rowIndex);
            this.data.CopyTo(clone.Data, 0);
            return clone;
        }

        #endregion
    }
}
