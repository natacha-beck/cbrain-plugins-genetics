
#
# CBRAIN Project
#
# Copyright (C) 2008-2023
# The Royal Institution for the Advancement of Learning
# McGill University
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# This module will extract pattern based on a regex
# specified in the descriptor by:
#
#   "cbrain:integrator_modules": {
#     "BoutiquesGroupFileByPattern": [ "input_id", "FileCollection", "A regex to extract prefix from the file name" ]
#   }
#
# It is possible to specify a FileCollection type in order to know on which
# input file to apply the regex.
#
# The extracted names will populate the input with value-choices
# according to what was extracted.
#
# Then one task will be created for each choice.
module BoutiquesGroupFileByPattern

    # Note: to access the revision info of the module,
    # you need to access the constant directly, the
    # object method revision_info() won't work.
    Revision_info=CbrainFileRevision[__FILE__] #:nodoc:

    # If more than one choice is detected in the FileCollection,
    # the input_id option will be populated with prefix
    # extracted based on pattern specified by the regex in the descriptor.
    def descriptor_for_form #:nodoc:
        descriptor             = super.dup()
        input_id, type, regex_string,  = descriptor.custom_module_info('BoutiquesGroupFileByPattern')

        # Get the userfile_ids from the params
        # Return immediately if there is no exactly one file.
        userfile_ids = params["interface_userfile_ids"] || []
        userfiles    =  Userfile.find(userfile_ids).select{|x| x.is_a?(type.constantize)}
        return descriptor if userfiles.size != 1

        # Extract the list of sample names from the file names
        userfile     = userfiles.first

        # Fill the input with the list of prefix
        regex      = Regexp.new(regex_string)
        file_names = userfile.provider_collection_index(:top, :regular).map(&:name)
        file_names.map! {|x| Pathname.new(x).basename.to_s }
        input      = descriptor.input_by_id(input_id)

        input["value-choices"] =  file_names.map do |f_n|
                                    f_n.match(regex) && Regexp.last_match[1]
                                  end.compact.uniq

        descriptor
    end

    # One task will be created by value of
    # the input specified in the descriptor.
    def final_task_list #:nodoc:
      descriptor     = self.descriptor_for_final_task_list
      input_id, _, _ = descriptor.custom_module_info('BoutiquesGroupFileByPattern')

      params_values = self.invoke_params[input_id]
      return super if params_values.blank? || params_values.size == 1

      # Create one task for each value
      params_values.map do |value|
        task = self.dup
        task.description    = task.description || ""
        task.description   += "\n\nRun with the following value #{value} for the input #{input_id}"
        # Add the value to the description... for information
        task.invoke_params[input_id] = [value]
        task
      end

    end
  end
