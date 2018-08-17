using System;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    /// <summary>
    /// The exception that is thrown when a linear algebra operation that requires symmetric positive definite matrices is 
    /// applied to an indefinite matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class AsymmetricMatrixException: IndefiniteMatrixException
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="AsymmetricMatrixException"/> class.
        /// </summary>
        public AsymmetricMatrixException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="AsymmetricMatrixException"/> class with a specified error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public AsymmetricMatrixException(string message) : base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="AsymmetricMatrixException"/> class with a specified error message 
        /// and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
        ///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
        public AsymmetricMatrixException(string message, Exception inner) : base(message, inner)
        { }
    }
}
