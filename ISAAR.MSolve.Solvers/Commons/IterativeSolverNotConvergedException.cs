using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Solvers.Commons
{
    /// <summary>
    /// The exception that is thrown when an iterative solver or a domain decomposition solver that uses iterative methods 
    /// cannot converge for a given problem and set of parameters. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class IterativeSolverNotConvergedException : Exception
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="IterativeSolverNotConvergedException"/> class.
        /// </summary>
        public IterativeSolverNotConvergedException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="IterativeSolverNotConvergedException"/> class with a specified 
        /// error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public IterativeSolverNotConvergedException(string message) : base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="IterativeSolverNotConvergedException"/> class with a specified 
        /// error message and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">
        /// The exception that is the cause of the current exception. If the innerException parameter is not a null reference, 
        /// the current exception is raised in a catch block that handles the inner exception. 
        /// </param>
        public IterativeSolverNotConvergedException(string message, Exception inner) : base(message, inner)
        { }
    }
}
