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
 * NormalDistribution.cs, 27.03.2007
 * 
 * 09.08.2006: Initial version
 * 17.08.2006: Renamed field my and property My to mu/Mu to consistently use english names
 * 21.09.2006: Adapted to change in base class (field "generator" declared private (formerly protected) 
 *               and made accessible through new protected property "Generator")
 * 06.03.2007: Defined the UpdateHelpers method, which is called in setters of the Mu and Sigma property
 *				 and ensures that the changed distribution parameter is taken into account when generating 
 *				 the next random number
 * 27.03.2007: Overridden the now virtual base class method Reset, so the NormalDistribution is properly reset
 *               in any case (previously an already computed random number may still be returned)
 * 
 */

#region original copyrights
//   -*- C++ -*-
/*****************************************************************************
 *
 *   |_|_|_  |_|_    |_    |_|_|_  |_		     C O M M U N I C A T I O N
 * |_        |_  |_  |_  |_        |_		               N E T W O R K S
 * |_        |_  |_  |_  |_        |_		                     C L A S S
 *   |_|_|_  |_    |_|_    |_|_|_  |_|_|_|_	                 L I B R A R Y
 *
 * $Id: Normal.c,v 1.2 2002/01/14 11:37:33 spee Exp $
 *
 * CNClass: CNNormal --- CNNormal (Gaussian) distributed random numbers
 *
 *****************************************************************************
 * Copyright (C) 1992-1996   Communication Networks
 *                           Aachen University of Technology
 *                           D-52056 Aachen
 *                           Germany
 *                           Email: cncl-adm@comnets.rwth-aachen.de
 *****************************************************************************
 * This file is part of the CN class library. All files marked with
 * this header are free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.  This library is
 * distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
 * License for more details.  You should have received a copy of the GNU
 * Library General Public License along with this library; if not, write
 * to the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139,
 * USA.
 *****************************************************************************
 * original Copyright:
 * -------------------
 * Copyright (C) 1988 Free Software Foundation
 *    written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 * 
 * This file is part of the GNU C++ Library.  This library is free
 * software; you can redistribute it and/or modify it under the terms of
 * the GNU Library General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your
 * option) any later version.  This library is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU Library General Public License for more details.
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *****************************************************************************/
#endregion

using System;
using Troschuetz.Random;

namespace Troschuetz.Random
{
	/// <summary>
	/// Provides generation of normal distributed random numbers.
	/// </summary>
	/// <remarks>
    /// The implementation of the <see cref="NormalDistribution"/> type bases upon information presented on
    ///   <a href="http://en.wikipedia.org/wiki/Normal_distribution">Wikipedia - Normal distribution</a>
    ///   and the implementation in the <a href="http://www.lkn.ei.tum.de/lehre/scn/cncl/doc/html/cncl_toc.html">
    ///   Communication Networks Class Library</a>.
    /// </remarks>
	public class NormalDistribution : Distribution
	{
		#region instance fields
        /// <summary>
        /// Gets or sets the parameter mu which is used for generation of normal distributed random numbers.
        /// </summary>
        /// <remarks>Call <see cref="IsValidMu"/> to determine whether a value is valid and therefor assignable.</remarks>
        public double Mu
        {
            get
            {
                return this.mu;
            }
            set
            {
                if (this.IsValidMu(value))
                {
                    this.mu = value;
					this.UpdateHelpers();
                }
            }
        }

        /// <summary>
        /// Stores the parameter mu which is used for generation of normal distributed random numbers.
        /// </summary>
        private double mu;

        /// <summary>
		/// Gets or sets the parameter sigma which is used for generation of normal distributed random numbers.
		/// </summary>
		/// <remarks>Call <see cref="IsValidSigma"/> to determine whether a value is valid and therefor assignable.</remarks>
		public double Sigma
		{
			get
			{
				return this.sigma;
			}
			set
			{
                if (this.IsValidSigma(value))
                {
                    this.sigma = value;
					this.UpdateHelpers();
                }
        	}
		}

		/// <summary>
		/// Stores the parameter sigma which is used for generation of normal distributed random numbers.
		/// </summary>
		private double sigma;

        /// <summary>
        /// Stores a precomputed normal distributed random number that will be returned the next time 
        ///   <see cref="NextDouble"/> gets called.
        /// </summary>
        /// <remarks>
        /// Two new normal distributed random numbers are generated every other call to <see cref="NextDouble"/>.
        /// </remarks>
        private double helper1;

        /// <summary>
        /// Stores a value indicating whether <see cref="NextDouble"/> was called twice since last generation of 
        ///   normal distributed random numbers.
        /// </summary>
        /// <remarks>
        /// Two new normal distributed random numbers are generated every other call to <see cref="NextDouble"/>.
        /// </remarks>
        private bool helper2;
		#endregion

		#region construction
		/// <summary>
        /// Initializes a new instance of the <see cref="NormalDistribution"/> class, using a 
        ///   <see cref="StandardGenerator"/> as underlying random number generator.
		/// </summary>
        public NormalDistribution()
            : this(new StandardGenerator())
		{
		}

		/// <summary>
        /// Initializes a new instance of the <see cref="NormalDistribution"/> class, using the specified 
        ///   <see cref="Generator"/> as underlying random number generator.
        /// </summary>
        /// <param name="generator">A <see cref="Generator"/> object.</param>
        /// <exception cref="ArgumentNullException">
        /// <paramref name="generator"/> is NULL (<see langword="Nothing"/> in Visual Basic).
        /// </exception>
        public NormalDistribution(Generator generator)
            : base(generator)
        {
            this.mu = 1.0;
            this.sigma = 1.0;

			this.UpdateHelpers();
        }
		#endregion
	
		#region instance methods
        /// <summary>
        /// Determines whether the specified value is valid for parameter <see cref="Mu"/>.
        /// </summary>
        /// <param name="value">The value to check.</param>
        /// <returns><see langword="true"/>.</returns>
        public bool IsValidMu(double value)
        {
            return true;
        }

        /// <summary>
        /// Determines whether the specified value is valid for parameter <see cref="Sigma"/>.
		/// </summary>
		/// <param name="value">The value to check.</param>
		/// <returns>
		/// <see langword="true"/> if value is greater than 0.0; otherwise, <see langword="false"/>.
		/// </returns>
		public bool IsValidSigma(double value)
		{
			return value > 0.0;
		}

		/// <summary>
		/// Updates the helper variables that store intermediate results for generation of normal distributed random 
		///   numbers.
		/// </summary>
		private void UpdateHelpers()
		{
			this.helper1 = 0.0;
			this.helper2 = false;
		}
		#endregion

		#region overridden Distribution members
		/// <summary>
		/// Resets the normal distribution, so that it produces the same random number sequence again.
		/// </summary>
		/// <returns>
		/// <see langword="true"/>, if the normal distribution was reset; otherwise, <see langword="false"/>.
		/// </returns>
		public override bool Reset()
		{
			bool result = base.Reset();
			if (result)
			{
				this.UpdateHelpers();
			}

			return result;
		}
		
		/// <summary>
		/// Gets the minimum possible value of normal distributed random numbers.
		/// </summary>
        public override double Minimum
		{
			get
			{
				return double.MinValue;
			}
		}

		/// <summary>
		/// Gets the maximum possible value of normal distributed random numbers.
		/// </summary>
        public override double Maximum
		{
			get
			{
				return double.MaxValue;
			}
		}

		/// <summary>
		/// Gets the mean value of normal distributed random numbers.
		/// </summary>
        public override double Mean
		{
			get
			{
				return this.mu;
			}
		}
		
		/// <summary>
		/// Gets the median of normal distributed random numbers.
		/// </summary>
        public override double Median
		{
			get
			{
                return this.mu;
			}
		}
		
		/// <summary>
		/// Gets the variance of normal distributed random numbers.
		/// </summary>
        public override double Variance
		{
			get
			{
				return Math.Pow(this.sigma, 2.0);
			}
		}
		
		/// <summary>
		/// Gets the mode of normal distributed random numbers.
		/// </summary>
        public override double[] Mode
		{
			get
			{
                return new double[] { this.mu };
			}
		}
		
		/// <summary>
		/// Returns a normal distributed floating point random number.
		/// </summary>
		/// <returns>A normal distributed double-precision floating point number.</returns>
        public override double NextDouble()
        {
            if (this.helper2)
            {
                this.helper2 = false;

                return this.helper1;
            }
            else
            {
				this.helper2 = true;
				
				while (true)
                {
                    double v1 = 2.0 * this.Generator.NextDouble() - 1.0;
                    double v2 = 2.0 * this.Generator.NextDouble() - 1.0;
                    double w = v1 * v1 + v2 * v2;

                    if (w <= 1)
                    {
                        double y = Math.Sqrt(-2.0 * Math.Log(w) / w) * this.sigma;
                        this.helper1 = v2 * y + this.mu;
                        return v1 * y + this.mu;
                    }
                }
            }
        }
		#endregion
	}
}