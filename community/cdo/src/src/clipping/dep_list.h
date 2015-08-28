/**
 * @file dep_list.h
 * @brief Utlity functions for dependency lists 
 *
 * Dependency lists and small general functions
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://redmine.dkrz.de/doc/YAC/html/index.html
 *
 * This file is part of YAC.
 *
 * YAC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with YAC.  If not, see <http://www.gnu.org/licenses/gpl.txt>.
 */

#ifndef DEP_LIST_H
#define DEP_LIST_H

/** \example test_dep_list.c
 * This contains examples on how to use struct dep_list.
 */

/** \brief dependency list
 *
 * data structure for storing and working with dependencies between unsigned integer data\n
 * possible applications of a dep_list are:
 *    - the corner ids for all cells (cell 0 -> corners 0, 1, 4, 3; cell 1 -> corners 1, 2, 5, 4; 2 -> 3, 4 , 7, 6; 3 -> 4, 5, 8, 7)
 *    - neighbour cell ids (cell 0 -> neighbours 1, 2, 3; cell 1 -> neighbours 0, 2, 3; 2 -> 0, 1, 3; 3 -> 0, 1, 2)
 *
 * The element index (in the example the cell index) goes from 0 to n-1, where n is the total number of elements in the dependency list.\n
 * The dependencies can have any valid unsigned integer number. However, one has to be careful when using \ref yac_invert_dep_list, because if the values of the dependencies are to big the resulting dependency list can be big.
 */
struct dep_list {
   unsigned num_elements;           //!< total number of elements in the dependency list
   unsigned * num_deps_per_element; //!< array containing the number of dependencies per element
   unsigned * dependencies;         //!< array containing all dependencies

   /**
    * \brief exclusive prefix sum for faster access of dependencies array
    *
    * this array is automatically generated and contains a exclusive prefix sum\n
    *   [0]=0, [1]=num_dep[0], [2]=num_dep[0]+num_dep[1], ...\n
    * see: http://en.wikipedia.org/wiki/Prefix_sum
    */
   unsigned * prescan;
};

/**
 * initialises a dependency list, can be used in order to avoid memory access violations
 * @param[in,out] list dependency list that is to be initialised
 */
void yac_init_dep_list (struct dep_list * list);

/**
 * initialises an empty dependency list, can be used if the dependencies themselves are not know a priori
 * @param[in,out] list dependency list that is to be initialised
 * @param[in] num_elements number of elements in the dependency list
 */
void yac_init_empty_dep_list(struct dep_list * list, unsigned num_elements);

/**
 * initialises dependency list and sets the dependencies (a previous call to init_dep_list is not mandatory)
 * @param[in,out] list dependency list that is to be initialised
 * @param[in] num_elements number of elements in the dependency list
 * @param[in] num_deps_per_element array of size num_elements that contains for each element the number of associated dependencies
 * @param[in] dependencies array that contains all dependencies, entries 0 to num_deps_per_element[0]-1 contain dependencies for element 0, entries num_deps_per_element[0] to num_deps_per_element[0]+num_deps_per_element[1]-1 contain the dependencies for element 1, etc.
 *
 * \remark This routine makes no copy of the array passed to it. The array are assumed to be on the heap and will be freed by a call to \ref yac_free_dep_list.
 */
void yac_set_dependencies (struct dep_list * list, unsigned num_elements,
                           unsigned * num_deps_per_element, unsigned * dependencies);
/**
 * adds dependencies to an existing dependency list
 * @param[in,out] list dependency list to be edited
 * @param[in] element element of the list to be changed
 * @param[in] num_dependencies number of dependencies to be added
 * @param[in] dependencies array containing the dependencies, which are supposed to be added to the list
 */
void yac_add_dependencies (struct dep_list * list, unsigned element,
                           unsigned num_dependencies, unsigned * dependencies);

/**
 * gets all dependencies for a given element
 * @param[in] list dependency list
 * @param[in] index element index
 * @return pointer to a constant array containing the dependencies for the given element index
 *
 * \remark the number of dependencies for the respective element can be retrieved from list.num_deps_per_element[index]
 */
unsigned const * yac_get_dependencies_of_element (struct dep_list list, unsigned index);

/**
 * gets the total number of dependencies stored in the given list
 * @param[in] list dependency list
 * @return total number of dependencies in list
 */
unsigned yac_get_total_num_dependencies(struct dep_list list);

/**
 * gets the position of a dependency for a given element in the dep_list::dependencies array
 * @param[in] list dependency list
 * @param[in] index element index
 * @param[in] dependency dependency of provided element
 * @return position of dependency in dep_list::dependencies array \n -1 in case index and/or dependency is invalid
 */
unsigned yac_get_dependency_index(struct dep_list list, unsigned index, unsigned dependency);

/**
 * gets the position of the first dependency of an element in the dep_list::dependencies array
 * @param[in] list dependency list
 * @param[in] index element index
 * @return position of the first dependency of an element in the dep_list::dependencies array
 */
unsigned yac_get_dependency_offset(struct dep_list list, unsigned index);

/**
 * search for a given dependency in a dependency list
 * @param[in] list dependency list
 * @param[in] dependency dependency that is to be searched for
 * @return 0 if the list does not contain the respective dependency
 */
unsigned yac_list_contains_dependency(struct dep_list list, unsigned dependency);

/**
 * gets the element and dependency associated to a given position in the dep_list::dependencies array
 * @param[in] list dependency list
 * @param[in] dep_index position in dep_list::dependencies
 * @param[out] index element index associated to dep_index
 * @param[out] dependency dependency associated to dep_index
 */
void yac_get_dependency(struct dep_list list, unsigned dep_index,
                        unsigned * index, unsigned * dependency);

/**
 * generates an inverted dependency list
 * @param[in] dep dependency list that is to be inverted
 * @param[out] inv_dep inverted version of dep (initialising inv_dep is no required)
 */
void yac_invert_dep_list(struct dep_list dep, struct dep_list * inv_dep);

/**
 * removes all dependencies of the provided elements
 * @param[in,out] dep dependency list
 * @param[in] element_indices array with all element indices for which the dependencies have to be removed from dep
 * @param[in] num_elements number of indices in element_indices array
 * \remark element indices == -1 are being ignored
 */
void yac_remove_dependencies_of_elements(struct dep_list * dep, unsigned * element_indices,
                                         unsigned num_elements);

/**
 * removes all given dependencies
 * @param[in,out] dep dependency list
 * @param[in] dependencies array containing the dependencies that have to be removed
 * @param[in] num_dependencies number of dependencies in dependencies array
 * \remark dependencies == -1 are being ignored
 */
void yac_remove_dependencies(struct dep_list * dep, unsigned * dependencies,
                             unsigned num_dependencies);

/**
 * makes a copy of a dependency list
 * @param[in] src source dependency list
 * @param[out] tgt target dependency list (initialising tgt is not required)
 */
void yac_copy_dep_list(struct dep_list src, struct dep_list * tgt);

/**
 * packs a dependency list into a buffer that can be sent to other processes via the MPI library for example
 * @param[in] list dependency list to be packed
 * @param[in,out] buf pointer to an array into which the dependency list is to packed (if too small this array will be reallocated)
 * @param[in] offset offset in the buffer at which the data is to be written
 * @param[out] data_size size of the data the is added to the buffer
 * @param[in,out] buf_size size of the buffer (this is being updated if necessary)
 * \remark *buf == NULL and *buf_size=0 are valid input values
 * @see \ref yac_unpack_dep_list
 */
void yac_pack_dep_list(struct dep_list list, unsigned ** buf, unsigned offset,
                       unsigned * data_size, unsigned * buf_size);
/**
 * unpacks a packed dependency list
 * @param[out] list dependency into which the unpacked data is to be written (initialising is not required)
 * @param[in] buf buffer containing packed dependency list
 * @param[out] data_size size of the data the was occupied by the packed dependency list
 * @see yac_pack_dep_list
 */
void yac_unpack_dep_list(struct dep_list * list, unsigned * buf, unsigned * data_size);

/**
 * frees a dependency list
 * @param[in,out] list dependency list to be freed
 */
void yac_free_dep_list(struct dep_list * list);

#endif // DEP_LIST_H
