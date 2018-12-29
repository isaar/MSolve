using System;

namespace ISAAR.MSolve.Solvers.Commons
{
    /// <summary>
    /// The exception that is thrown when the chosen solver does not match the discrete model. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InvalidSolverException: Exception 
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidSolverException"/> class.
        /// </summary>
        public InvalidSolverException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidSolverException"/> class with a specified error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public InvalidSolverException(string message) : base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidSolverException"/> class with a specified error message 
        /// and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
        ///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
        public InvalidSolverException(string message, Exception inner) : base(message, inner)
        { }
    }
}
