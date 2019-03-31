/** \file md5.h
 *
 * \brief Contains functions to calculate MD5 sums
 */
#pragma once

/// calculate md5 sum of a data block in memory
std::string md5 (void *data, int number_of_bytes);
/// calculate md5 sum of a file on disk
std::string md5 (const std::string &filename);
