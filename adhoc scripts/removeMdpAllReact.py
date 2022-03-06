# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 15:37:02 2022

@author: new
"""
# work in progress


def removeMdpAllReact(self, vReact, vTar):
    """Removes the tar data of a reaction from the RDML data.

    Args:
        self: The class self parameter.
        vReact: The reaction id.
        vTar: The target id.

    Returns:
        Nothing, updates RDML data.
    """

    reacts = _get_all_children(self._node, "react")
    for react in reacts:
        react_datas = _get_all_children(react, "data")
        for react_data in react_datas:
            mdpdat = _get_first_child(react_data, "mdp")
        #     if forId is not None:
        #         if forId.attrib['id'] == vTar:
        #             react.remove(react_data)
        #             break
        # remain_data += len(_get_all_children(react, "data"))
        # partit = _get_first_child(react, "partitions")
        # if partit is not None:
        #     remain_data += len(_get_all_children(partit, "data"))
        # if remain_data == 0:
        #     self._node.remove(react)
        #     break
    return