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
 * ChiSquareDistribution.cs, 27.03.2007
 * 
 * 17.08.2006: Initial version
 * 27.03.2007: Overridden the now virtual base class method Reset, so the ChiSquareDistribution is properly reset
 *               in any case (must explicitely reset the underlying NormalDistribution)
 * 
 */

using System;
using Troschuetz.Random;

namespace Troschuetz.Random
{
	/// <summary>
    /// Provides generation of chi-square distributed random numbers.
	/// </summary>
	/// <remarks>
	/// The implementation of the <see cref="ChiSquareDistribution"/> type bases upon information presented on
    ///   <a href="http://en.wikipedia.org/wiki/Chi-square_distribution">Wikipedia - Chi-square distribution</a>.
    /// </remarks>
	public class ChiSquareDistribution : Distribution
	{
		#region instance fields
		/// <summary>
        /// Gets or sets the parameter alpha which is used for generation of chi-square distributed random numbers.
		/// </summary>
		/// <remarks>Call <see cref="IsValidAlpha"/> to determine whether a value is valid and therefor assignable.</remarks>
		public int Alpha
		{
			get
			{
                return this.alpha;
			}
			set
			{
                if (this.IsValidAlpha(value))
                {
                    this.alpha = value;
                }
        	}
		}

		/// <summary>
        /// Stores the parameter alpha which is used for generation of chi-square distributed random numbers.
		/// </summary>
        private int alpha;

        /// <summary>
        /// Stores a <see cref="NormalDistribution"/> object used for generation of chi-square distributed random numbers.
        /// </summary>
        private NormalDistribution normalDistribution;
        #endregion

		#region construction
		/// <summary>
        /// Initializes a new instance of the <see cref="ChiSquareDistribution"/> class, using a 
        ///   <see cref="StandardGenerator"/> as underlying random number generator.
		/// </summary>
        public ChiSquareDistribution()
            : this(new StandardGenerator())
		{
		}
		
		/// <summary>
        /// Initializes a new instance of the <see cref="ChiSquareDistribution"/> class, using the specified 
        ///   <see cref="Generator"/> as underlying random number generator.
        /// </summary>
        /// <param name="generator">A <see cref="Generator"/> object.</param>
        /// <exception cref="ArgumentNullException">
        /// <paramref name="generator"/> is NULL (<see langword="Nothing"/> in Visual Basic).
        /// </exception>
        public ChiSquareDistribution(Generator generator)
            : base(generator)
        {
            this.alpha = 1;
            this.normalDistribution = new NormalDistribution(generator);
            this.normalDistribution.Mu = 0.0;
            this.normalDistribution.Sigma = 1.0;
        }
		#endregion
	
		#region instance methods
		/// <summary>
        /// Determines whether the specified value is valid for parameter <see cref="Alpha"/>.
		/// </summary>
		/// <param name="value">The value to check.</param>
		/// <returns>
		/// <see langword="true"/> if value is greater than 0; otherwise, <see langword="false"/>.
		/// </returns>
        public bool IsValidAlpha(int value)
		{
			return value > 0;
		}
        #endregion

		#region overridden Distribution members
		/// <summary>
		/// Resets the chi-square distribution, so that it produces the same random number sequence again.
		/// </summary>
		/// <returns>
		/// <see langword="true"/>, if the chi-square distribution was reset; otherwise, <see langword="false"/>.
		/// </returns>
		public override bool Reset()
		{
			bool result = base.Reset();
			if (result)
			{
				result = this.normalDistribution.Reset();
			}

			return result;
		}

		/// <summary>
        /// Gets the minimum possible value of chi-square distributed random numbers.
		/// </summary>
        public override double Minimum
		{
			get
			{
				return 0.0;
			}
		}

		/// <summary>
        /// Gets the maximum possible value of chi-square distributed random numbers.
		/// </summary>
        public override double Maximum
		{
			get
			{
				return double.MaxValue;
			}
		}

		/// <summary>
        /// Gets the mean value of chi-square distributed random numbers.
		/// </summary>
        public override double Mean
		{
			get
			{
                return this.alpha;
			}
		}
		
		/// <summary>
        /// Gets the median of chi-square distributed random numbers.
		/// </summary>
        public override double Median
		{
			get
			{
				return this.alpha - 2.0 / 3.0;
			}
		}
		
		/// <summary>
        /// Gets the variance of chi-square distributed random numbers.
		/// </summary>
        public override double Variance
		{
			get
			{
                return 2.0 * this.alpha;
			}
		}
		
		/// <summary>
        /// Gets the mode of chi-square distributed random numbers.
		/// </summary>
        public override double[] Mode
		{
            get
            {
                if (this.alpha >= 2)
                {
                    return new double[] { this.alpha - 2.0 };
                }
                else
                {
                    return new double[] { };
                }
            }
		}
		
		/// <summary>
        /// Returns a chi-square distributed floating point random number.
		/// </summary>
        /// <returns>A chi-square distributed double-precision floating point number.</returns>
        public override double NextDouble()
		{
            double sum = 0.0;
            for (int i = 0; i < this.alpha; i++)
            {
                sum += Math.Pow(this.normalDistribution.NextDouble(), 2);
            }

            return sum;
		}
        #endregion
    }
}