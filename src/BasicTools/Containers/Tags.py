# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

from __future__ import annotations # python 3 compatibility

from typing import Iterable, Optional,  List, Tuple, Iterator, Collection

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject, froze_it
from BasicTools.NumpyDefs import PBasicIndexType, ArrayLike

@froze_it
class Tag(BaseOutputObject):
    """A Tag is an object to store a name and ids.
    internals it has a buffer to easily fill the
    content incrementally.
    """
    def __init__(self,tagname):
        super().__init__()

        self.name = tagname
        self._id = np.empty(0,dtype=PBasicIndexType)
        self.cpt = 0

    def __eq__(self, other:Tag) -> bool:
        """Equal operator, return True only if names and ids are equal,
        internally this function calls Tighten

        Parameters
        ----------
        other : object
            the other tag

        Returns
        -------
        bool
            True if names and GetIds of self and other are equal
            False otherwise
        """
        if not isinstance(other, Tag):
            return False

        if self.name != other.name:
            return False

        self.Tighten()
        other.Tighten()

        if not np.array_equal(self._id,other._id):
            return False

        return True

    def AddToTag(self,tid: int | Collection):
        """Add id or ids to the list of tags
        more memory is allocated if needed

        Parameters
        ----------
        tid : PBasicIndexType or Iterable
            the id or ids to added to the tag.
        """
        if isinstance(tid, Collection):
            if len(self._id) <= self.cpt+len(tid):
                self._id = np.resize(self._id, (len(self._id)*2+len(tid),))
            self._id[self.cpt:self.cpt+len(tid)] = tid
            self.cpt += len(tid)
        else:
            if len(self._id) <= self.cpt:
                self._id = np.resize(self._id, (self.cpt*2+1,))

            self._id[self.cpt] = tid
            self.cpt += 1

    def __len__(self) -> int :
        """return the number of ids in this tag. The internal memory
        can be larger

        Returns
        -------
        int
            number of ids
        """
        return self.cpt

    def Tighten(self)->None:
        """Release non used memory
        """
        if self._id.shape[0] != self.cpt:
            self._id = np.resize(self._id, (self.cpt,))

    def RemoveDoubles(self)->None:
        """Remove doubles and release non used memory
        """
        self.SetIds(self._id[0:self.cpt])

    def SetIds(self, ids:ArrayLike):
        """Set the ids of this tag, a copy is made and a
        RemoveDoubles is executed to ensure ids are present only once.
        old data is lost.

        Parameters
        ----------
        ids : ArrayLike
            the list of ids
        """
        self._id = np.unique(np.asarray(ids,dtype=PBasicIndexType))
        self.cpt = len(self._id)

    def SetId(self, pos:int, idd:int):
        """Set the value of the id in the position pos.
        The user is responsible to allocate the memory first

        Parameters
        ----------
        pos : int
            position
        idd : int
            the value
        """
        self._id[pos] = idd

    def Allocate(self, allocationSize:int):
        """Allocate the memory for n objects
        the user is responsible of filling each individual position.
        Internally numpy.resize is used.

        Parameters
        ----------
        l : int
            the number of element to allocate
        """
        self.cpt = allocationSize
        self.Tighten()

    def GetIds(self) -> np.ndarray:
        """Return the Ids in the tag

        Returns
        -------
        np.ndarray
            _description_
        """
        self.Tighten()
        return self._id

    def GetIdsAsMask(self,totalNumberOfObjects:int=0, output:np.ndarray=None, erase:bool=True)-> np.ndarray:
        """Generate a numpy array of dtype=bool of size totalNumberOfObjects with the indexes of this
        as True

        if output is not None, this array will be uses as output.

        Parameters
        ----------
        totalNumberOfObjects : int, optional
            total number of Objects to allocate the output
            This argument is not compatible with the output argument
            by default None
        output : bool, optional
            output array for work in place
            This argument is not compatible with the totalNumberOfObjects
            by default None
        erase : bool, optional
            set all values to False before working on output,
            by default True

        Returns
        -------
        np.ndarray
            Array of bool
        """

        self.Tighten()

        if output is None:
            output = np.zeros(totalNumberOfObjects, dtype=bool)
        else:
            if erase :
                output.fill(False)

        output[self._id] = True

        return output

    def Merge(self, other:Tag):
        """Merge the other tag into this

        Parameters
        ----------
        other : Tag, optional
            _description_, by default None
        """
        if other is not None:
            self.AddToTag(other.GetIds())

        self.RemoveDoubles()

    def __str__(self) -> str:
        res = ''
        res  = str(self.name) + " " + str(self.cpt) + " "
        return res

    def Copy(self) -> Tag:
        """Create a copy of this Tag

        Returns
        -------
        Tag
            a Tag with the same name and a copy of the ids
        """
        tag = Tag(self.name)
        tag._id = np.copy(self.GetIds())
        tag.cpt = len(tag._id)
        return tag

class Tags(BaseOutputObject):
    """Storage to hold tag :
        It has a Dict like interface
        Propagate some function to every tag
        Functions to create and delete tag
    """
    def __init__(self):
        super(Tags,self).__init__()
        self.storage = []

    def Tighten(self)->None:
        """Call Tag.Tighten on every tag """

        for tag in self:
            tag.Tighten()

    def RemoveDoubles(self)->None:
        """Call Tag.RemoveDoubles on every tag """
        for tag in self:
            tag.RemoveDoubles()

    def AddTag(self, item:Tag ) -> Tag:
        """Add a tag to the container

        Parameters
        ----------
        item : Tag
            The tag to be added, this tag must contain a unique name

        Returns
        -------
        Tag
            the Tag added

        Raises
        ------
        Exception
            if a tag with the same name exist already
        """

        if item.name in self:
            raise Exception("Cant add the tag two times!!")# pragma: no cover

        self.storage.append(item)
        return item

    def DeleteTags(self, tagNames:List[str]):
        """Ensure that tags with the names in the tagNames are deleted

        Parameters
        ----------
        tagNames : List[str]
            List of tagnames
        """

        self.storage = [ tag for tag in self.storage if tag.name not in tagNames ]

    def CreateTag(self, name: str, errorIfAlreadyCreated:Optional[bool]=True) -> Tag:
        """Create a new tag with the name "name" and return it

        Parameters
        ----------
        name : str
            Name of the tag to be created
        errorIfAlreadyCreated : Optional[bool], optional
            If False and a tag with the same name exit, then this tag is returned , by default True

        Returns
        -------
        Tag
            Tag with the name : name

        Raises
        ------
        Exception
            If errorIfAlreadyCreated is True and the tag already exist, a Exception is raised
        """
        if name in self:
            if errorIfAlreadyCreated :
                raise Exception("Tag name '"+name+"' already exist")# pragma: no cover
            else:
                return self[name]
        else:
            return self.AddTag(Tag(name))

    def RenameTag(self, name: str, newName: str, noError:Optional[bool]= False):
        """Rename a tag if exist, nonError controls if a exception is raised

        Parameters
        ----------
        name : str
            Name of the tag to change the name
        newName : str
            New name
        noError : Optional[bool], optional
            Option to silently do nothing if no tag exist with name "name" , by default False

        Raises
        ------
        Exception
            if noError is False (default) a Exception is raised if a tag with the name "name" does not exist
        """
        if name in self:
            self[name].name = newName
        else:
            if noError:
                return
            raise Exception("Tag '" + str(name) + "' does not exist")# pragma: no cover

    def RemoveEmptyTags(self):
        """Remove empty tag
        """
        for tagname in  list(self.keys()):
            tag = self[tagname]
            if tag.cpt == 0:
                self.storage.remove(tag)

## function to act like a dict

    def keys(self):
        """Dict API, the keys are the names"""
        return [ x.name for x in self.storage ]

    def __contains__(self, tagName=str):
        """Dict API, check if a tag with the name tagName exist"""
        for item in self.storage:
            if item.name == tagName:
                return True
        return False

    def __getitem__(self, tagName:str) -> Tag:
        """Dict API, Return the tag with the tagName

        Parameters
        ----------
        tagName : str
            Tag name to extract

        Returns
        -------
        Tag
            Tag whit the name tagName

        Raises
        ------
        KeyError
            if no tag with the name tagName exist
        """

        for item in self.storage:
            if item.name == tagName:
                return item
        raise KeyError("Tag '"+ str(tagName) + "' not found")# pragma: no cover

    def __iter__(self) -> Iterator[Tag]:
        """Dict API, iterable

        Returns
        -------
        Iterable[Tag]
            Tags inside this container
        """
        return iter(self.storage)

    def __len__(self) -> int:
        """Return the number of tags

        Returns
        -------
        int
            number of tags
        """

        return len(self.storage)

    def items(self) -> Iterable[Tuple[str,Tag]]:
        """Dict API

        Returns
        -------
        Iterable[Tuple[str,Tag]]
            _description_
        """
        return [(v.name,v) for v in self.storage]

    def __str__(self) -> str:
        res = ''
        res  = str(list(self.keys()))
        return res

    def Copy(self)-> Tags:
        """Make a new copy of the tag containers.

        Returns
        -------
        Tags
            A copy containing copies of each tag.
        """
        res = Tags()
        for tag in self:
            res.AddTag(tag.Copy())
        return res

def CheckIntegrity():

    tag = Tag("Tag")
    tag.Allocate(2)
    tag.SetId(0,1)
    tag.SetId(1,2)

    if len(tag) != 2: # pragma: no cover
        raise Exception("Incorrect number of element in the container")

    if np.any(tag.GetIds() !=  [1,2]):# pragma: no cover
        raise Exception("Error in the tag")

    tagII = Tag("Tag")
    tagII.AddToTag([1])
    tagII.AddToTag(2)
    tagII.AddToTag(2)
    tagII.RemoveDoubles()

    if tag != tagII:# pragma: no cover
        print(tag)
        print(tagII)
        raise Exception("tags are not equal")

    tagII.name = "tagII"
    if tag == tagII:# pragma: no cover
        raise Exception("equal but with different names")

    if tag == Tag('Tag'):# pragma: no cover
        raise Exception("equal but with different ids")

    tagII.SetIds([1,2])

    if np.any(tagII.GetIdsAsMask(3,) != [False, True,True]):# pragma: no cover
        raise Exception("Error in the generation of the mask")

    out= np.random.random(3)

    if np.any(tagII.GetIdsAsMask(3, output=out, erase=True ) != [False, True,True]):# pragma: no cover
        raise Exception("Error in the generation of the mask")

    tag.Merge(tagII)

    if len(tag) != 2:# pragma: no cover
        raise Exception("Error in the merge")

    print(tag)
    print(tag.Copy())

    tags = Tags()
    tags.AddTag(tag)
    tags.AddTag(Tag("tag2"))
    tags.AddTag(Tag("tag3"))
    tags.Tighten()

    if len(tags) != 3: # pragma: no cover
        raise Exception("Error in the number of tags ")

    tags.DeleteTags(["newtoto","tag"])

    tagsII = Tags()
    tagsII.AddTag(Tag("tag2"))
    tagsII.CreateTag("tag3")
    tagsII.CreateTag("tag3", errorIfAlreadyCreated=False)

    try:
        tagsII.CreateTag("tag3", errorIfAlreadyCreated=True)
    except Exception:
        pass
    else: # pragma: no cover
        raise Exception("This function must fail")


    tagsII.RenameTag("tag2","TAG2")
    tagsII.RenameTag("tag2","TAG2", noError=True)

    tagsII.RemoveDoubles()
    tagsII.RemoveEmptyTags()
    if len(tagsII) !=0: # pragma: no cover
        raise Exception("Unable to remove the empty tags")

    print(tags.Copy())
    # Dict interface
    print(tags.items())


    return "OK"

if __name__ == '__main__': # pragma: no cover
    print(CheckIntegrity())
