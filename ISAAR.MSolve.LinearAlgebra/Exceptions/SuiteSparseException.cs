using System;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    /// <summary>
    /// The exception that is thrown when a call to SuiteSparse library fails.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SuiteSparseException : Exception
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="SuiteSparseException"/> class.
        /// </summary>
        public SuiteSparseException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="SuiteSparseException"/> class with a specified error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public SuiteSparseException(string message) : base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="SuiteSparseException"/> class with a specified error message 
        /// and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
        ///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
        public SuiteSparseException(string message, Exception inner) : base(message, inner)
        { }
    }
}
