/* Copyright (C) 2025 Jonathan Parkinson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CDR_ASSIGNMENT_UTILITIES_HEADER_H
#define CDR_ASSIGNMENT_UTILITIES_HEADER_H

// C++ headers
#include <vector>
#include <tuple>
#include <string>

// Library headers

// Project headers
#include "../annotator_classes/numbering_constants.inc"


namespace CDRConversionUtilities {


// We have to be able to convert between CDRs as defined in
// one scheme and the corresponding numbering in another. This
// requires a table of "magic numbers" that define the start
// and end of CDRs in different schemes. Unfortunately this
// is lengthy because of the need to convert between each
// pair of schemes. These magic numbers appear below.

// ************** NUMBERED WITH IMGT ************************ //
static constexpr std::array<int, 6> IMGT_CDR_BREAKPOINTS =
{27, 39, 56, 66, 105, 118};
static constexpr std::array<int, 6> IMGT_NMBR_KABAT_H_CDR_BREAKPOINTS =
{32, 41, 55, 75, 107, 118};
static constexpr std::array<int, 6> IMGT_NMBR_KABAT_L_CDR_BREAKPOINTS =
{24, 41, 56, 70, 105, 118};
static constexpr std::array<int, 6> IMGT_NMBR_MARTIN_H_CDR_BREAKPOINTS =
{27, 38, 57, 65, 107, 118};
static constexpr std::array<int, 6> IMGT_NMBR_MARTIN_L_CDR_BREAKPOINTS =
{26, 39, 56, 66, 107, 117};
static constexpr std::array<int, 6> IMGT_NMBR_AHO_CDR_BREAKPOINTS =
{25, 39, 56, 76, 107, 117};
static constexpr std::array<int, 6> IMGT_NMBR_NORTH_H_CDR_BREAKPOINTS =
{24, 41, 55, 67, 105, 118};
static constexpr std::array<int, 6> IMGT_NMBR_NORTH_L_CDR_BREAKPOINTS =
{24, 41, 55, 70, 105, 118};

// ************** NUMBERED WITH AHO ************************ //
static constexpr std::array<int, 6> AHO_CDR_BREAKPOINTS =
{25, 41, 58, 78, 109, 138};
static constexpr std::array<int, 6> AHO_NMBR_KABAT_H_CDR_BREAKPOINTS =
{32, 43, 57, 77, 109, 139};
static constexpr std::array<int, 6> AHO_NMBR_KABAT_L_CDR_BREAKPOINTS =
{24, 43, 58, 72, 107, 139};
static constexpr std::array<int, 6> AHO_NMBR_MARTIN_H_CDR_BREAKPOINTS =
{27, 40, 59, 67, 109, 139};
static constexpr std::array<int, 6> AHO_NMBR_MARTIN_L_CDR_BREAKPOINTS =
{26, 41, 58, 68, 109, 138};
static constexpr std::array<int, 6> AHO_NMBR_IMGT_CDR_BREAKPOINTS =
{27, 41, 58, 68, 107, 139};
static constexpr std::array<int, 6> AHO_NMBR_NORTH_H_CDR_BREAKPOINTS =
{24, 43, 57, 69, 107, 139};
static constexpr std::array<int, 6> AHO_NMBR_NORTH_L_CDR_BREAKPOINTS =
{24, 43, 57, 72, 107, 139};


// ************** NUMBERED WITH KABAT ************************ //
static constexpr std::array<int, 6> KABAT_HEAVY_CDR_BREAKPOINTS =
{31, 36, 50, 66, 95, 103};
static constexpr std::array<int, 6> KABAT_LIGHT_CDR_BREAKPOINTS =
{24, 35, 50, 57, 89, 98};
static constexpr std::array<int, 6> KABAT_NMBR_MARTIN_H_CDR_BREAKPOINTS =
{26, 36, 52, 57, 95, 103};
static constexpr std::array<int, 6> KABAT_NMBR_MARTIN_L_CDR_BREAKPOINTS =
{26, 33, 50, 53, 91, 97};
static constexpr std::array<int, 6> KABAT_NMBR_IMGT_H_CDR_BREAKPOINTS =
{26, 36, 51, 58, 93, 103};
static constexpr std::array<int, 6> KABAT_NMBR_IMGT_L_CDR_BREAKPOINTS =
{27, 33, 50, 53, 89, 98};
static constexpr std::array<int, 6> KABAT_NMBR_AHO_H_CDR_BREAKPOINTS =
{24, 36, 51, 67, 95, 102};
static constexpr std::array<int, 6> KABAT_NMBR_AHO_L_CDR_BREAKPOINTS =
{25, 33, 50, 62, 91, 97};
static constexpr std::array<int, 6> KABAT_NMBR_NORTH_H_CDR_BREAKPOINTS =
{23, 36, 50, 59, 93, 103};
static constexpr std::array<int, 6> KABAT_NMBR_NORTH_L_CDR_BREAKPOINTS =
{24, 35, 49, 57, 89, 98};




// ************** NUMBERED WITH MARTIN ************************ //
static constexpr std::array<int, 6> MARTIN_HEAVY_CDR_BREAKPOINTS =
{26, 33, 52, 57, 95, 103};
static constexpr std::array<int, 6> MARTIN_LIGHT_CDR_BREAKPOINTS =
{26, 33, 50, 53, 91, 97};
static constexpr std::array<int, 6> MARTIN_NMBR_KABAT_H_CDR_BREAKPOINTS =
{31, 36, 50, 66, 95, 103};
static constexpr std::array<int, 6> MARTIN_NMBR_KABAT_L_CDR_BREAKPOINTS =
{24, 35, 50, 57, 89, 98};
static constexpr std::array<int, 6> MARTIN_NMBR_IMGT_H_CDR_BREAKPOINTS =
{26, 34, 51, 58, 93, 103};
static constexpr std::array<int, 6> MARTIN_NMBR_IMGT_L_CDR_BREAKPOINTS =
{27, 33, 50, 53, 89, 98};
static constexpr std::array<int, 6> MARTIN_NMBR_AHO_H_CDR_BREAKPOINTS =
{24, 34, 51, 67, 95, 102};
static constexpr std::array<int, 6> MARTIN_NMBR_AHO_L_CDR_BREAKPOINTS =
{25, 33, 50, 62, 91, 97};
static constexpr std::array<int, 6> MARTIN_NMBR_NORTH_H_CDR_BREAKPOINTS =
{23, 36, 50, 59, 93, 103};
static constexpr std::array<int, 6> MARTIN_NMBR_NORTH_L_CDR_BREAKPOINTS =
{24, 35, 49, 57, 89, 98};



/// @brief Assigns cdr labels to an input sequence.
/// @param numbering A list of numbering codes valid for the
///        scheme specified when the class was created. This is
///        the first element of the tuple returned by 'analyze_seq'.
/// @param chain A valid chain ('H', 'K', 'L'). 'K' and 'L' are treated
///        as equivalent.
/// @param cdr_labels A vector in which the output will be stored.
/// @param scheme The scheme -- one of 'imgt', 'aho', 'martin', 'kabat' --
///        that was used to number the antibodies.
/// @param cdr_scheme The scheme -- one of 'imgt', 'aho', 'martin', 'kabat',
///        'north' -- that should be used to construct CDRs. This may
///        or may not be the same as scheme.
void assign_cdr_labels(const std::vector<std::string> &numbering,
        const std::string &chain,
        std::vector<std::string> &cdr_labeling,
        const std::string &scheme,
        const std::string &cdr_scheme);

}  // namespace SequenceUtilities


#endif
