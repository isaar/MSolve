using System;

namespace ISAAR.MSolve.Solvers.Commons
{
    /// <summary>
    /// The exception that is thrown when a solver or global matrix assembler does not match the analysis model or the format of 
    /// the linear system matrix. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InvalidMatrixFormatException: Exception
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidMatrixFormatException"/> class.
        /// </summary>
        public InvalidMatrixFormatException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidMatrixFormatException"/> class with a specified error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public InvalidMatrixFormatException(string message) : base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidMatrixFormatException"/> class with a specified error message 
        /// and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
        ///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
        public InvalidMatrixFormatException(string message, Exception inner) : base(message, inner)
        { }
    }
}
