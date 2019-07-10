using System;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    /// <summary>
    /// The exception that is thrown when a linear algebra operation that requires invertible matrices is applied to a singular
    /// matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SingularMatrixException: IndefiniteMatrixException
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="SingularMatrixException"/> class.
        /// </summary>
        public SingularMatrixException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="SingularMatrixException"/> class with a specified error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public SingularMatrixException(string message): base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="SingularMatrixException"/> class with a specified error message 
        /// and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
        ///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
        public SingularMatrixException(string message, Exception inner): base(message, inner)
        { }
    }
}
