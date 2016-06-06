/*
 * Copyright © 2006 Stefan Troschütz (stefan@troschuetz.de)
 * 
 * This file is part of Troschuetz.Random Class Library.
 * 
 * Troschuetz.Random is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * 
 * StudentsTDistribution.cs, 27.03.2007
 * 
 * 17.08.2006: Initial version
 * 27.03.2007: Overridden the now virtual base class method Reset, so the StudentsTDistribution is properly reset
 *               in any case (must explicitely reset the underlying NormalDistribution)
 * 
 */

using System;
using Troschuetz.Random;

namespace Troschuetz.Random
{
	/// <summary>
    /// Provides generation of t-distributed random numbers.
	/// </summary>
	/// <remarks>
	/// The implementation of the <see cref="StudentsTDistribution"/> type bases upon information presented on
    ///   <a href="http://en.wikipedia.org/wiki/Student%27s_t-distribution">Wikipedia - Student's t-distribution</a> and
    ///   <a href="http://www.xycoon.com/stt_random.htm">Xycoon - Student t Distribution</a>.
    /// </remarks>
	public class StudentsTDistribution : Distribution
	{
		#region instance fields
		/// <summary>
        /// Gets or sets the parameter nu which is used for generation of t-distributed random numbers.
		/// </summary>
        /// <remarks>Call <see cref="IsValidNu"/> to determine whether a value is valid and therefor assignable.</remarks>
        public int Nu
		{
			get
			{
                return this.nu;
			}
			set
			{
                if (this.IsValidNu(value))
                {
                    this.nu = value;
                    this.UpdateHelpers();
                }
        	}
		}

		/// <summary>
        /// Stores the parameter nu which is used for generation of t-distributed random numbers.
		/// </summary>
        private int nu;

        /// <summary>
        /// Stores a <see cref="NormalDistribution"/> object used for generation of t-distributed random numbers.
        /// </summary>
        private NormalDistribution normalDistribution;

        /// <summary>
        /// Stores a <see cref="ChiSquareDistribution"/> object used for generation of t-distributed random numbers.
        /// </summary>
        private ChiSquareDistribution chiSquareDistribution;
        #endregion

		#region construction
		/// <summary>
        /// Initializes a new instance of the <see cref="StudentsTDistribution"/> class, using a 
        ///   <see cref="StandardGenerator"/> as underlying random number generator.
		/// </summary>
        public StudentsTDistribution()
            : this(new StandardGenerator())
		{
		}
		
		/// <summary>
        /// Initializes a new instance of the <see cref="StudentsTDistribution"/> class, using the specified 
        ///   <see cref="Generator"/> as underlying random number generator.
        /// </summary>
        /// <param name="generator">A <see cref="Generator"/> object.</param>
        /// <exception cref="ArgumentNullException">
        /// <paramref name="generator"/> is NULL (<see langword="Nothing"/> in Visual Basic).
        /// </exception>
        public StudentsTDistribution(Generator generator)
            : base(generator)
        {
            this.nu = 1;
            this.normalDistribution = new NormalDistribution(generator);
            this.normalDistribution.Mu = 0.0;
            this.normalDistribution.Sigma = 1.0;
            this.chiSquareDistribution = new ChiSquareDistribution(generator);
            this.UpdateHelpers();
        }
		#endregion
	
		#region instance methods
		/// <summary>
        /// Determines whether the specified value is valid for parameter <see cref="Nu"/>.
		/// </summary>
		/// <param name="value">The value to check.</param>
		/// <returns>
		/// <see langword="true"/> if value is greater than 0; otherwise, <see langword="false"/>.
		/// </returns>
        public bool IsValidNu(int value)
		{
			return value > 0;
		}

        /// <summary>
        /// Updates the helper variables that store intermediate results for generation of t-distributed random 
        ///   numbers.
        /// </summary>
        private void UpdateHelpers()
        {
            this.chiSquareDistribution.Alpha = this.nu;
        }
        #endregion

		#region overridden Distribution members
		/// <summary>
		/// Resets the Student's t-distribution, so that it produces the same random number sequence again.
		/// </summary>
		/// <returns>
		/// <see langword="true"/>, if the Student's t-distribution was reset; otherwise, <see langword="false"/>.
		/// </returns>
		public override bool Reset()
		{
			bool result = base.Reset();
			if (result)
			{
				result = this.normalDistribution.Reset();
				if (result)
				{
					result = this.chiSquareDistribution.Reset();
				}
			}

			return result;
		}

		/// <summary>
        /// Gets the minimum possible value of t-distributed random numbers.
		/// </summary>
        public override double Minimum
		{
			get
			{
				return double.MinValue;
			}
		}

		/// <summary>
        /// Gets the maximum possible value of t-distributed random numbers.
		/// </summary>
        public override double Maximum
		{
			get
			{
				return double.MaxValue;
			}
		}

        /// <summary>
        /// Gets the mean value of t-distributed random numbers.
		/// </summary>
        public override double Mean
		{
			get
			{
                if (this.nu > 1)
                {
                    return 0.0;
                }
                else
                {
                    return double.NaN;
                }
			}
		}
		
		/// <summary>
        /// Gets the median of t-distributed random numbers.
		/// </summary>
        public override double Median
		{
			get
			{
				return 0.0;
			}
		}
		
		/// <summary>
        /// Gets the variance of t-distributed random numbers.
		/// </summary>
        public override double Variance
		{
			get
			{
                if (this.nu > 2)
                {
                    return (double)this.nu / ((double)this.nu - 2.0);
                }
                else
                {
                    return double.NaN;
                }
			}
		}
		
		/// <summary>
        /// Gets the mode of t-distributed random numbers.
		/// </summary>
        public override double[] Mode
		{
            get
            {
                return new double[] { 0.0 };
            }
		}
		
		/// <summary>
        /// Returns a t-distributed floating point random number.
		/// </summary>
        /// <returns>A t-distributed double-precision floating point number.</returns>
        public override double NextDouble()
		{
            return this.normalDistribution.NextDouble() / Math.Sqrt(this.chiSquareDistribution.NextDouble() / (double)this.nu);
		}
        #endregion
    }
}