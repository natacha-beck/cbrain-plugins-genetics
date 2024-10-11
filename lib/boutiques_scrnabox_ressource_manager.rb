
#
# CBRAIN Project
#
# Copyright (C) 2008-2022
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

module BoutiquesScrnaboxRessourceManager

  # Note: to access the revision info of the module,
  # you need to access the constant directly, the
  # object method revision_info() won't work.
  Revision_info=CbrainFileRevision[__FILE__] #:nodoc:

  def descriptor_for_cluster_commands
    descriptor                = super.dup

    suggested_resources_by_steps = descriptor.custom_module_info("BoutiquesScrnaboxRessourceManager")

    invoke_steps          = invoke_params["steps"] # '1' or '2' or '1-8'
    first_step, last_step = invoke_steps.split("-").map(&:to_i)
    last_step           ||= first_step
    selected_steps        = (first_step..last_step).to_a.map(&:to_s)

    walltime = 0
    ram      = 0

    # Use suggested resources for each step
    suggested_resources_by_steps.each do |step_option, resources|
        step, option = step_option.split("_")

        next unless selected_steps.include?(step)
        next if option.present? && invoke_params[option].blank?

        step_wall = resources["walltime"]
        walltime += step_wall

        step_ram  = resources["ram"]
        ram       = step_ram if step_ram > ram
    end



    cb_error "Cannot compute time or ram limits from selected steps" if ram == 0 || walltime == 0

    descriptor["suggested-resources"]["walltime-estimate"] = walltime * 60
    descriptor["suggested-resources"]["ram"]               = ram

    return descriptor
  end

end

