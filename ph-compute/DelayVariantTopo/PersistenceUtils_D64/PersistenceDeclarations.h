#pragma once

#include "PersistenceUtils_D64.h"
#include "PersistenceBarcodes.h"

typedef float FType;

typedef boost::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<FType>> CPersistenceBarcodesPtr;
typedef boost::shared_ptr<const NPersistenceUtils::CPersistenceBarcodes<FType>> CPersistenceBarcodesCPtr;

typedef std::vector<CPersistenceBarcodesPtr> TypeBarcodesPtrVec;