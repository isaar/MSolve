using System;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    /// <summary>
    /// The exception that is thrown when a linear algebra operation cannot be executed due to the arguments (vectors or 
    /// matrices) having incompatible dimension (number of rows or columns).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NonMatchingDimensionsException: Exception
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="NonMatchingDimensionsException"/> class.
        /// </summary>
        public NonMatchingDimensionsException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="NonMatchingDimensionsException"/> class with a specified error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public NonMatchingDimensionsException(string message): base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="NonMatchingDimensionsException"/> class with a specified error message 
        /// and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
        ///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
        public NonMatchingDimensionsException(string message, Exception inner): base(message, inner)
        { }
    }
}
