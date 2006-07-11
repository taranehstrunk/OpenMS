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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PERSISTENTOBJECT_H
#define OPENMS_FORMAT_PERSISTENTOBJECT_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

namespace OpenMS
{ 
	class PersistenceManager;

  /**
  	@brief Base class for all persistent objects.
  	
  	Interface for all classes that can be stored persistently e.g. in DB or text file.
  	The storage itself is handled by a PersistenceManager.
  	
  	@ingroup Format
  */
  class PersistentObject
  {

    public:
      /// Default constructor
      PersistentObject();

      /// Destructor
      virtual ~PersistentObject();
			
      /// Assignment operator
      PersistentObject& operator= (const PersistentObject& rhs);

      /**
      	@brief Method for writing an object to a stream
      	
      	
      */
      virtual void persistentWrite(PersistenceManager& pm, const char* name=0) const throw (Exception::Base) =0;

      /**
      	@brief Method for reading an object from a stream
      	
      	
      */
      virtual void persistentRead(PersistenceManager& pm) throw (Exception::Base) =0;
			
      /**
      	@brief Returns the persistence id
      	
      	This id is only used by some kinds of PersistenceManagers.
      	E.g. in the DBAdapter the id is used to connect the object to the data stored in the DB.
      */
      const UID& getPersistenceId() const;

      /// Returns a reference the persistence id
      UID& getPersistenceId();

      /**
      	@brief Sets the persistence id
      	
      	This id is only used by some kinds of PersistenceManagers.
      	E.g. in the DBAdapter the id is used to connect the object to the data stored in the DB.
      	<BR>
      	Do not set the persistence id unless you know what you are doing!
      */
      void setPersistenceId(const UID& persistence_id);

      /**
      	@brief Clears the persistence id
      	
      	Sets the id to 0.<br>
      	@param deep determines which ids are cleared. <tt>false</tt> means that only the id of the current object is reset. 
      	<tt>true</tt> means that the ids of all sub-objects are reset as well (default).
      */
      void clearId(bool deep = true);

    protected:
    
      ///A persistence id used to refer the data back to the source
      UID persistence_id_;

      /**
      	@brief Clears the persistence id of all sub-objects.
      	
      	
      */
      virtual void clearChildIds_() =0;
  };

}
#endif
