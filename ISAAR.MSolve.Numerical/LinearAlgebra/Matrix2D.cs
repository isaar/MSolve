using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class Matrix2D : IMatrix2D, ILinearlyCombinable
    {
        private bool isTransposed = false;
        private int rows, columns;
        private double[,] data;

        public Matrix2D(int rows, int columns)
        {
            this.rows = rows;
            this.columns = columns;
            data = new double[rows,columns];
        }

        public Matrix2D(double[,] data)
        {
            this.rows = data.GetLength(0);
            this.columns = data.GetLength(1);
            this.data = data;
        }

        public double[,] Data
        {
            get { return data; }
        }

        #region IMatrix2D Members

        public int Rows
        {
            get { return isTransposed ? columns : rows; }
        }

        public int Columns
        {
            get { return isTransposed ? rows : columns; }
        }

        public double this[int x, int y]
        {
            get 
            {
                return isTransposed ? data[y, x] : data[x, y]; 
            }
            set 
            {
                if (isTransposed)
                    data[y, x] = value;
                else
                    data[x, y] = value; 
            }
        }

        /// <summary>
        /// this = this + scalar * matrix
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="scalar"></param>
        public void AxpyIntoThis(IMatrix2D matrix, double scalar)
        {
            if ((this.Rows != matrix.Rows) || (this.Columns != matrix.Columns))
            {
                throw new ArgumentException("The two matrices must have the same dimensions.");
            }
            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < Columns; ++j)
                {
                    this[i, j] += scalar * matrix[i, j];
                }
            }
        }

        public void Multiply(IVector vIn, double[] vOut)
        {
            if (isTransposed == false)
                MultiplyNormal(vIn, vOut);
            else
                MultiplyTranspose(vIn, vOut);
        }

        public void Multiply(double[] vIn, double[] vOut)
        {
            Matrix2D AA = new Matrix2D(data);
            for (int i = 0; i < rows; i++)
            {
                vOut[i] = 0;
                for (int j = 0; j < columns; j++)
                    vOut[i] += AA.Data[i, j] * vIn[j];
            }
        }

        private void MultiplyNormal(IVector vIn, double[] vOut)
        {
            Matrix2D AA = new Matrix2D(data);
            for (int i = 0; i < rows; i++)
            {
                vOut[i] = 0;
                for (int j = 0; j < columns; j++)
                    vOut[i] += AA.Data[i, j] * vIn[j];
            }
        }

        private void MultiplyTranspose(IVector vIn, double[] vOut)
        {
            Matrix2D AA = new Matrix2D(data);
            for (int j = 0; j < columns; j++)
            {
                vOut[j] = 0;
                for (int i = 0; i < rows; i++)
                    vOut[j] += AA.Data[i, j] * vIn[i]; 
            }
        }

        public void Scale(double scale)
        {
            double[,] mData = data;
            for (int i = 0; i < mData.GetLength(0); i++)
                for (int j = 0; j < mData.GetLength(1); j++)
                    mData[i, j] *= scale;
        }

        #endregion

        public void Clear()
        {
            Array.Clear(data, 0, data.Length);
        }

        public Matrix2D Transpose()
        {
            var m = new Matrix2D(this.data);
            m.isTransposed = true;
            return m;
        }

        public static Matrix2D operator *(Matrix2D A, Matrix2D B)
        {
            if (A.Columns != B.Rows) throw new ArgumentException("Matrix sizes mismatch.");

            double[,] c = new double[A.Rows, B.Columns];
            Matrix2D AA = new Matrix2D(A.data);
            AA.isTransposed = A.isTransposed;

            for (int i = 0; i < A.Rows; i++)
                for (int k = 0; k < B.Columns; k++)
                    for (int j = 0; j < B.Rows; j++)
                        c[i, k] += AA[i, j] * B[j, k];

            //var sw = System.IO.File.CreateText(@"d:\BeamTransformedSecondPart.txt");
            //for (int i = 0; i < A.Rows; i++)
            //{
            //    var s = string.Empty;
            //    for (int j = 0; j < B.Columns; j++)
            //        s += c[i, j].ToString() + ";";
            //    sw.WriteLine(s);
            //}
            //sw.Close();
            return new Matrix2D(c);
        }

		public void Add(Matrix2D B)
		{
			
			double[,] mData = data;
			if (mData.GetLength(0) != B.rows&& mData.GetLength(1) != B.Columns) throw new ArgumentException("Matrix sizes mismatch.");
			for (int i = 0; i < mData.GetLength(0); i++)
			for (int j = 0; j < mData.GetLength(1); j++)
				mData[i, j] += B[i,j];
		}

		public static Vector operator *(Matrix2D A, Vector b)
        {
            if (A.Columns != b.Length) throw new ArgumentException("Matrix sizes mismatch.");

            double[] c = new double[A.Rows];
            Matrix2D AA = new Matrix2D(A.data);
            AA.isTransposed = A.isTransposed;

            for (int i = 0; i < A.Rows; i++)
                for (int j = 0; j < A.Columns; j++)
                    c[i] += AA[i, j] * b[j];
            return new Vector(c);
        }

        public static Matrix2D operator *(Matrix2D A, double b)
        {
            Matrix2D AA = new Matrix2D(A.data as double[,]);
            for (int i = 0; i < A.Rows; i++)
                for (int j = 0; j < A.Columns; j++)
                    AA[i, j] = AA[i, j] * b;
            return AA;
        }

        public void SVD(double[] w, double[,] v)
        {
            //      double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
            //      double precision dsqrt,dmax1,dabs,dsign

            int i, j, k, l, ii, i1, kk, k1, ll, l1, mn, its;
            double c, f, g, h, s, x, y, z, scale, anorm, myTemp;
            bool bEsc1;

            int m = rows;
            int n = columns;
            double[,] u = new double[rows, columns];
            double[,] a = data as double[,];
            bool matu = false;
            bool matv = true;
            double[] rv1 = new double[n];

            //  ierr = 0;

            Array.Copy(a, u, a.GetLength(0) * a.GetLength(1)); 
            //a.CopyTo(u, 0);
            //for (i = 1; i <= m; i++)
            //    for (j = 1; j <= n; j++)
            //        u[i - 1, j - 1] = a[i - 1, j - 1];

            //     .......... householder reduction to bidiagonal form ..........
            g = 0;
            scale = 0;
            anorm = 0;
            l = 0;
            for (i = 1; i <= n; i++)
            {
                l = i + 1;
                rv1[i - 1] = scale * g;
                g = 0;
                s = 0;
                scale = 0;

                if (i <= m)
                {
                    for (k = i; k <= m; k++)
                        scale = scale + Math.Abs(u[k - 1, i - 1]);

                    if (scale != 0)
                    {
                        for (k = i; k <= m; k++)
                        {
                            u[k - 1, i - 1] = u[k - 1, i - 1] / scale;
                            s = s + Math.Pow(u[k - 1, i - 1], 2);
                        }

                        f = u[i - 1, i - 1];
                        g = -(Math.Sqrt(s));
                        if (f < 0)
                            g = g * (-1);

                        h = f * g - s;
                        u[i - 1, i - 1] = f - g;

                        if (i != n)
                        {
                            for (j = l; j <= n; j++)
                            {
                                s = 0;
                                for (k = i; k <= m; k++)
                                    s = s + u[k - 1, i - 1] * u[k - 1, j - 1];
                                f = s / h;
                                for (k = i; k <= m; k++)
                                    u[k - 1, j - 1] = u[k - 1, j - 1] + f * u[k - 1, i - 1];
                            }
                        }

                        for (k = i; k <= m; k++)
                            u[k - 1, i - 1] = scale * u[k - 1, i - 1];
                    }
                }

                w[i - 1] = scale * g;
                g = 0;
                s = 0;
                scale = 0;

                if (!((i > m) || (i == n)))
                {
                    for (k = l; k <= n; k++)
                        scale = scale + Math.Abs(u[i - 1, k - 1]);

                    if (scale != 0)
                    {
                        for (k = l; k <= n; k++)
                        {
                            u[i - 1, k - 1] = u[i - 1, k - 1] / scale;
                            s = s + Math.Pow(u[i - 1, k - 1], 2);
                        }

                        f = u[i - 1, l - 1];
                        g = -(Math.Sqrt(s));
                        if (f < 0)
                            g = g * (-1);
                        h = f * g - s;
                        u[i - 1, l - 1] = f - g;

                        for (k = l; k <= n; k++)
                            rv1[k - 1] = u[i - 1, k - 1] / h;

                        if (i != m)
                        {
                            for (j = l; j <= m; j++)
                            {
                                s = 0;
                                for (k = l; k <= n; k++)
                                    s = s + u[j - 1, k - 1] * u[i - 1, k - 1];
                                for (k = l; k <= n; k++)
                                    u[j - 1, k - 1] = u[j - 1, k - 1] + s * rv1[k - 1];
                            }
                        }

                        for (k = l; k <= n; k++)
                            u[i - 1, k - 1] = scale * u[i - 1, k - 1];
                    }
                }

                myTemp = Math.Abs(w[i - 1]) + Math.Abs(rv1[i - 1]);
                if (anorm < myTemp)
                    anorm = myTemp;
            }

            //     .......... accumulation of right-hand transformations ..........
            if (matv)
            {
                //     .......... for i=n step -1 until 1 do -- ..........
                for (ii = 1; ii <= n; ii++)
                {
                    i = n + 1 - ii;

                    if (i != n)
                    {
                        if (g != 0)
                        {
                            for (j = l; j <= n; j++)
                                //     .......... double division avoids possible underflow ..........
                                v[j - 1, i - 1] = (u[i - 1, j - 1] / u[i - 1, l - 1]) / g;

                            for (j = l; j <= n; j++)
                            {
                                s = 0;
                                for (k = l; k <= n; k++)
                                    s = s + u[i - 1, k - 1] * v[k - 1, j - 1];

                                for (k = l; k <= n; k++)
                                    v[k - 1, j - 1] = v[k - 1, j - 1] + s * v[k - 1, i - 1];
                            }
                        }

                        for (j = l; j <= n; j++)
                        {
                            v[i - 1, j - 1] = 0;
                            v[j - 1, i - 1] = 0;
                        }
                    }

                    v[i - 1, i - 1] = 1;
                    g = rv1[i - 1];
                    l = i;
                }
            }

            //     .......... accumulation of left-hand transformations ..........
            if (matu)
            {
                //     ..........for i=min(m,n) step -1 until 1 do -- ..........
                mn = n;
                if (m < n) mn = m;

                for (ii = 1; ii <= mn; ii++)
                {
                    i = mn + 1 - ii;
                    l = i + 1;
                    g = w[i - 1];

                    if (i != n)
                        for (j = l; j <= n; j++)
                            u[i - 1, j - 1] = 0;

                    if (g != 0)
                    {
                        if (i != mn)
                        {
                            for (j = l; j <= n; j++)
                            {
                                s = 0;

                                for (k = l; k <= m; k++)
                                    s = s + u[k - 1, i - 1] * u[k - 1, j - 1];
                                //     .......... double division avoids possible underflow ..........
                                f = (s / u[i - 1, i - 1]) / g;

                                for (k = i; k <= m; k++)
                                    u[k - 1, j - 1] = u[k - 1, j - 1] + f * u[k - 1, i - 1];
                            }
                        }

                        for (j = i; j <= m; j++)
                            u[j - 1, i - 1] = u[j - 1, i - 1] / g;
                    }
                    else
                    {
                        for (j = i; j <= m; j++)
                            u[j - 1, i - 1] = 0;
                    }

                    u[i - 1, i - 1] = u[i - 1, i - 1] + 1;
                }
            }

            //     .......... diagonalization of the bidiagonal form ..........
            //     .......... for k=n step -1 until 1 do -- ..........
            for (kk = 1; kk <= n; kk++)
            {
                k1 = n - kk;
                k = k1 + 1;
                its = 0;

                //     .......... test for convergence ..........
                // 520    .......... test for splitting.
                //                for l=k step -1 until 1 do -- ..........
                //    while (l != k)
                l1 = 0;
                for (; ; )
                {
                    bEsc1 = false;
                    for (ll = 1; ll <= k; ll++)
                    {
                        l1 = k - ll;
                        l = l1 + 1;

                        if ((Math.Abs(rv1[l - 1]) + anorm) == anorm)
                        {
                            bEsc1 = true;
                            break;
                        }
                        //     .......... rv1(1) is always zero, so there is no exit
                        //                through the bottom of the loop ..........

                        if ((Math.Abs(w[l1 - 1]) + anorm) == anorm)
                            break;
                    }

                    if (!(bEsc1))
                    {
                        //     .......... cancellation of rv1(l) if l greater than 1 ..........
                        c = 0;
                        s = 1.0;

                        for (i = l; i <= k; i++)
                        {
                            f = s * rv1[i - 1];
                            rv1[i - 1] = c * rv1[i - 1];
                            if ((Math.Abs(f) + anorm) == anorm)
                                break;

                            g = w[i - 1];
                            h = Math.Sqrt(f * f + g * g);
                            w[i - 1] = h;
                            c = g / h;
                            s = -f / h;
                            if (matu)
                                for (j = 1; j <= m; j++)
                                {
                                    y = u[j - 1, l1 - 1];
                                    z = u[j - 1, i - 1];
                                    u[j - 1, l1 - 1] = y * c + z * s;
                                    u[j - 1, i - 1] = -y * s + z * c;
                                }
                        }
                    }

                    z = w[k - 1];
                    if (l == k)
                        break;
                    //     .......... shift from bottom 2 by 2 minor ..........
                    //      if (its .eq. 30) go to 1000
                    its = its + 1;
                    x = w[l - 1];
                    y = w[k1 - 1];
                    g = rv1[k1 - 1];
                    h = rv1[k - 1];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
                    g = Math.Sqrt(f * f + 1);
                    /*      myTemp = fabs(g);
                          if (f < 0)
                            myTemp = myTemp * (-1); */

                    if (f < 0)
                        myTemp = -g;
                    else
                        myTemp = g;
                    f = ((x - z) * (x + z) + h * (y / (f + myTemp) - h)) / x;
                    //     .......... next qr transformation ..........
                    c = 1;
                    s = 1;

                    for (i1 = l; i1 <= k1; i1++)
                    {
                        i = i1 + 1;
                        g = rv1[i - 1];
                        y = w[i - 1];
                        h = s * g;
                        g = c * g;
                        z = Math.Sqrt(f * f + h * h);
                        rv1[i1 - 1] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = -x * s + g * c;
                        h = y * s;
                        y = y * c;
                        if (matv)
                        {
                            for (j = 1; j <= n; j++)
                            {
                                x = v[j - 1, i1 - 1];
                                z = v[j - 1, i - 1];
                                v[j - 1, i1 - 1] = x * c + z * s;
                                v[j - 1, i - 1] = -x * s + z * c;
                            }
                        }

                        z = Math.Sqrt(f * f + h * h);
                        w[i1 - 1] = z;
                        //     .......... rotation can be arbitrary if z is zero ..........
                        if (z != 0)
                        {
                            c = f / z;
                            s = h / z;
                        }
                        f = c * g + s * y;
                        x = -s * g + c * y;
                        if (matu)
                        {
                            for (j = 1; j <= m; j++)
                            {
                                y = u[j - 1, i1 - 1];
                                z = u[j - 1, i - 1];
                                u[j - 1, i1 - 1] = y * c + z * s;
                                u[j - 1, i - 1] = -y * s + z * c;
                            }
                        }
                    }

                    rv1[l - 1] = 0;
                    rv1[k - 1] = f;
                    w[k - 1] = x;
                    //      go to 520
                }  // end while

                //     .......... convergence ..........
                if (z < 0)
                {
                    //     .......... w(k) is made non-negative ..........
                    w[k - 1] = -z;
                    if (matv)
                    {
                        for (j = 1; j <= n; j++)
                            v[j - 1, k - 1] = -v[j - 1, k - 1];
                    }
                }
            }


            //// Free up u
            //  if (!matu)
            //  {
            //    for (i = 0; i < rows; i++) delete u[i];
            //    delete u;
            //  } 
        }        

        public void LinearCombinationGOAT(IList<double> coefficients, IList<IMatrix2D> matrices)
        {
            if (coefficients.Count != matrices.Count)
                throw new ArgumentException("Coefficient and matrix sizes do not match");
            foreach (var m in matrices)
                if (m.Rows != this.Rows || m.Columns != this.Columns)
                    throw new ArgumentException("Matrix sizes do not match");
            var localData = data as double[,];
            var localCoefficients = coefficients as IList<double>;

            for (int k = 0; k < coefficients.Count; k++)
            {
                var m = (Matrix2D)matrices[k];
                var d = localCoefficients[k];
                for (int i = 0; i < this.Rows; i++)
                    for (int j = 0; j < this.Columns; j++)
                        localData[i, j] += d * m.Data[i, j];
            }
        }

        public void LinearCombination(IList<double> coefficients, IList<IMatrix2D> matrices)
        {
            if (coefficients.Count != matrices.Count)
                throw new ArgumentException(String.Format("Coefficients and matrices count mismatch ({0} <> {1}).", coefficients.Count, matrices.Count));
            for (int i = 0; i < matrices.Count; i++)
                if (matrices[i].Rows != rows)
                    throw new ArgumentException(String.Format("Matrix at pos {0} has {1} rows instead of {2}.", i, matrices[i].Rows, rows));

            var cs = (IList<double>)coefficients;
            double[,] newData = new double[rows, columns];
            for (int k = 0; k < matrices.Count; k++)
            {
                var m = matrices[k];
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < columns; j++)
                        newData[i, j] += cs[k] * m[i, j];
                //for (int j = 0; j < data.Length; j++)
                //    newData[j] += cs[i] * m.Data[j];
            }

            Array.Copy(newData, data, data.Length);
        }

        public static Matrix2D FromVector(double[] vector)
        {
            var m = new Matrix2D(vector.Length, 1);
            for (int i = 0; i < m.Data.Length; i++)
                m.Data[i, 0] = vector[i];
            return m;
        }

        public static Matrix2D FromVectorTranspose(double[] vector)
        {
            var m = new Matrix2D(1, vector.Length);
            for (int i = 0; i < m.Data.Length; i++)
                m.Data[0, i] = vector[i];
            return m;
        }        
    }
}
