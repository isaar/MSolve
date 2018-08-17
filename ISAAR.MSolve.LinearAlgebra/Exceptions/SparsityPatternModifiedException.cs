using System;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    /// <summary>
    /// The exception that is thrown when modifying a zero entry that is not explicitly stored in the vector or matrix storage 
    /// format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SparsityPatternModifiedException: PatternModifiedException
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="SparsityPatternModifiedException"/> class.
        /// </summary>
        public SparsityPatternModifiedException()
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="SparsityPatternModifiedException"/> class with a specified error message.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        public SparsityPatternModifiedException(string message) : base(message)
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="SparsityPatternModifiedException"/> class with a specified error message 
        /// and a reference to the inner exception that is the cause of this exception.
        /// </summary>
        /// <param name="message">The error message that explains the reason for the exception.</param>
        /// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
        ///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
        public SparsityPatternModifiedException(string message, Exception inner) : base(message, inner)
        { }
    }
}
