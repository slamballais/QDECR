// The shape of a docs sidebar, shared by the guide and the reference.

export interface SidebarItem {
  href: string
  label: string
  /** Set the label as code: function names in the reference. */
  code?: boolean
}

export interface SidebarGroup {
  label: string
  /** A page for the group as a whole, which its label then links to. */
  href?: string
  /** List the items as numbered steps. */
  numbered?: boolean
  items: SidebarItem[]
}
