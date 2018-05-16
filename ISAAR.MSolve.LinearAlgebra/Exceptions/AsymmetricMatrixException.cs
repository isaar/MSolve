using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    public class AsymmetricMatrixException: IndefiniteMatrixException
    {
        public AsymmetricMatrixException()
        { }

        public AsymmetricMatrixException(string message) : base(message)
        { }

        public AsymmetricMatrixException(string message, Exception inner) : base(message, inner)
        { }
    }
}
