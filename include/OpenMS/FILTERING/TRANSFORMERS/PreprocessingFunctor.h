// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H
#define OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H

#include <vector>
#include <map>
#include <string>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{
	/**
    @defgroup Filtering Filtering
    
    @defgroup SpectraPreprocessing Spectra Preprocessors

    @ingroup Filtering
  */

  /**
  	@brief Base class for MSSpectrum preprocessing classes

		@ingroup SpectraPreprocessing
  */
  class PreprocessingFunctor : public FactoryProduct
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    PreprocessingFunctor();
		
    /// copy constructor
    PreprocessingFunctor(const PreprocessingFunctor& source);
		
    /// destructor
    virtual ~PreprocessingFunctor() {}
		// @}
		
		// @name Operators
		// @{
    /// assignment operator
    PreprocessingFunctor& operator=(const PreprocessingFunctor& source);
		// @}

		// @name Accessors
		// @{
		/// interface definition of the functor classes
		template <typename SpectrumType> void apply(SpectrumType& /* spectrum */) {}
		// @}
	};

}
#endif // OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H
